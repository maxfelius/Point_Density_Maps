'''
@author: Max Felius
@date: 30/09/2020

Script to create a resolution map. It will combine datasets and count the amount of points per pre-defined cell

TODO:
    - Create a better coordinates transformation scheme
    - Create saves in between steps
    - Create github repo for version control (mess right now) 

ISSUE:
    - Not all the datasets have pnt_rdx and pnt_rdy columns
'''

#imports
import numpy as np
import os
import sys
import time
import pandas as pd
from scipy import spatial
import progressbar

#import custom class
import rijksdriehoek

#main class
class create_resolution_map:
    def __init__(self,filename,resolution,filename_inter_data='combined_points.csv',save_intermediate_results=True):
        '''
        Initiate class

        filename should be a string or a list. filename should be the absolute path of the file

        :type filename: string/list[]
        :type resolution: int
        :rtype: None
        '''

        #Object variables
        self.filename = filename
        self.header = ['pnt_lon','pnt_lat','pnt_rdx','pnt_rdy']
        self.data = None
        self.xrange = None
        self.yrange = None
        self.xv = None
        self.yv = None
        self.radius = None
        self.kdtree = None
        self.cell_radius = resolution/2
        self.grid_counter = None
        self.grid_polygon = []
        self.save_intermediate_results = save_intermediate_results
        self.save_data = filename_inter_data

    def start(self):
        '''
        method to start the processing steps
        '''
        #1. load data
        self.data = self.Load_Data() #only columns for the rd_x and rd_y are returned

        #2. create evaluation grid
        self.Create_Grid()
        
        #3. Create KD-Tree
        self.kdtree = self.Create_KDTree()

        #4. Compute number of points per grid
        self.Count_Points()

        #5. Save the results to a good output
        return self.Create_df()
    
    def Load_Data(self):
        '''
        Method for loading the data
        '''
        #check if the intermediate file already exists
        if not os.path.exists(os.path.join('intermediate_data',self.save_data)):  

            if isinstance(self.filename,list):
                data = pd.DataFrame(columns=self.header)
                
                #for file_in in self.filename:
                #TEST: progressbar
                print('\nStarted reading the data...')
                for i in progressbar.progressbar(range(len(self.filename))):
                    file_in = self.filename[i]
                    temp = pd.read_csv(file_in)

                    if 'pnt_rdx' not in list(temp) or 'pnt_rdy' not in list(temp):
                        #creating the rd coordinates using lat/lon as input
                        rd = rijksdriehoek.Rijksdriehoek()
                        rd.from_wgs(temp['pnt_lat'],temp['pnt_lon'])

                        #create new dataframe swith rd coordinates
                        df = pd.DataFrame(np.array([rd.rd_x,rd.rd_y]).T,columns=['pnt_rdx','pnt_rdy'])
                        temp = temp.join(df)

                        data = data.append(temp[self.header], ignore_index=True)                
                    else:
                        data = data.append(temp[self.header], ignore_index=True)
                    print(f'\nAppended {file_in} to the dataset')
            else:
                data = pd.read_csv(self.filename)

            # save intermediate results such that the datasets don't have to be combined every new run
            if self.save_intermediate_results:
                if not os.path.isdir('intermediate_data'):
                    os.mkdir('intermediate_data')

                data[self.header].to_csv(os.path.join('intermediate_data',self.save_data))
            return data
        
        else:
            print(f'Reading pre-computed data {self.save_data}')

            return pd.read_csv(os.path.join('intermediate_data',self.save_data))

    def Create_Grid(self):
        '''
        Step 2 -> create grid
        '''
        max_x = max(self.data['pnt_rdx'])
        min_x = min(self.data['pnt_rdx'])

        max_y = max(self.data['pnt_rdy'])
        min_y = min(self.data['pnt_rdy'])

        resolution = self.cell_radius

        self.xrange = np.arange(min_x+resolution,max_x,resolution*2)
        self.yrange = np.arange(min_y+resolution,max_y,resolution*2)

        self.xv, self.yv = np.meshgrid(self.xrange,self.yrange)

    def Create_KDTree(self):
        '''
        Init the kdtree
        '''
        return spatial.cKDTree(self.data[['pnt_rdx','pnt_rdy']].values) #input is an array
    
    def Count_Points(self):
        '''
        Compute the grid cell

        TODO:
        '''
        #time keeping
        start = time.time()

        #A cell is a square. From the center to a corner the radius should be multiplied by sqrt(2).
        #The search radius is made 1% larger just in case.
        radius = np.sqrt(2)*self.cell_radius*1.01

        grid_counter = np.zeros((len(self.yrange),len(self.xrange)))

        print('Started counting points per cell...')    
        idx_y = 0
        idx_x = 0
        for idx in progressbar.progressbar(range(len(self.yrange)*len(self.xrange))):
            
            #update idx_y and yrange if needed
            if idx == idx_y*len(self.xrange):
                y = self.yrange[idx_y]
                idx_y += 1

            #update idx_x and xrange if needed
            if idx == 0:
                pass
            else:
                idx_x += 1

            if idx % len(self.xrange) == 0:
                idx_x = 0

            x = self.xrange[idx_x]

            #step 1
            subset = self.kdtree.query_ball_point([x,y],r=radius)
            # subset = [] #testing if the polygon is created correctly
            #step 2

            if not subset:
                pass
            else:
                values = np.array(self.data[['pnt_rdx','pnt_rdy']].iloc[subset].values)

                rdx = values[:,0]
                rdy = values[:,1]

                #created a function so it can potentially be accelerated by using a GPU (not implemented)
                def counter(rdx_in,rdy_in,cell_radius,x,y):
                    count = 0
                    for rdx,rdy in zip(rdx_in,rdy_in):
                        if x-cell_radius < rdx and x+cell_radius > rdx and y-cell_radius < rdy and y+cell_radius > rdy:
                            count += 1
                    return count

                count = counter(rdx,rdy,self.cell_radius,x,y)                

                grid_counter[idx_y-1,idx_x] = count
            
            #create the square polygon
            lbrd = rijksdriehoek.Rijksdriehoek(x-self.cell_radius,y-self.cell_radius)
            rbrd = rijksdriehoek.Rijksdriehoek(x+self.cell_radius,y-self.cell_radius)
            rtrd = rijksdriehoek.Rijksdriehoek(x+self.cell_radius,y+self.cell_radius)
            ltrd = rijksdriehoek.Rijksdriehoek(x-self.cell_radius,y+self.cell_radius)

            leftbot = lbrd.to_wgs()
            rightbot = rbrd.to_wgs()
            righttop = rtrd.to_wgs()
            lefttop = ltrd.to_wgs()

            self.grid_polygon.append(f'POLYGON (({leftbot[1]} {leftbot[0]},{rightbot[1]} {rightbot[0]},{righttop[1]} {righttop[0]},{lefttop[1]} {lefttop[0]}))')

        self.grid_counter = grid_counter
        print(f'Elapsed time is {time.time()-start} seconds...')

    def Create_df(self):
        '''
        Output a dataframe that can be saved
        '''
        X = np.ravel(self.xv)
        Y = np.ravel(self.yv)
        Z = np.ravel(self.grid_counter)

        rd = rijksdriehoek.Rijksdriehoek(X,Y)

        # wgs = np.array(list(map(rijksdriehoek.rd_to_wgs,X,Y)))
        wgs = np.array(rd.to_wgs()).T

        lat = wgs[:,0]
        lon = wgs[:,1]

        data_out = np.array([X,Y,lon,lat,Z,self.grid_polygon])

        return pd.DataFrame(data_out.T,columns=['rd_x','rd_y','lon','lat','Counts','wkt'])    
