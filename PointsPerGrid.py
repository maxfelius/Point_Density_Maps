'''
@author: Max Felius

Script build for the sinkhole workshop on 24 and 25th of October.

Version 1

'''

#imports
import numpy as np
import matplotlib.pyplot as plt
from decouple import config
import os
import pandas as pd
from scipy import spatial
import time
import rijksdriehoek

class pointCounter:
    def __init__(self,filename):
        '''
        Initialize the class
        
        TODO:
            1. Load Data (DONE)
            2. Create evaluation grid
            3. Create KD-Tree Object
            4. Compute number of points per grid point
            5. Save results to a pandas df
            6. Create notebook for plotting the results

        '''
        #Object variables
        self.filename = filename
        self.data = None
        self.xrange = None
        self.yrange = None
        self.xv = None
        self.yv = None
        self.radius = None
        self.kdtree = None
        self.resolution = 100 #[m] - Actual resolution is twice this (i.e. resolution = radius)
        self.grid_counter = None
        self.grid_polygon = []

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
        

    def Create_Points(self):
        '''
        method to test the error in the conversion
        '''
        X = self.data['pnt_rdx']
        Y = self.data['pnt_rdy']
        
        rd = rijksdriehoek.Rijksdriehoek(X,Y)

        # wgs = np.array(list(map(rijksdriehoek.rd_to_wgs,X,Y)))
        wgs = np.array(rd.to_wgs()).T
    
        lat = wgs[:,0]
        lon = wgs[:,1]

        coordinates = np.array([lon,lat])

        return pd.DataFrame(coordinates.T,columns=['lon','lat'])

    def Load_Data(self):
        '''
        Step 1 -> loading data
        '''
        path = config('main')
        filename = self.filename
        data = pd.read_csv(os.path.join(path,filename))

        return data[['pnt_rdx','pnt_rdy']]
    
    def Create_Grid(self):
        '''
        Step 2 -> create grid
        '''
        max_x = max(self.data['pnt_rdx'])
        min_x = min(self.data['pnt_rdx'])

        max_y = max(self.data['pnt_rdy'])
        min_y = min(self.data['pnt_rdy'])

        resolution = self.resolution

        self.xrange = np.arange(min_x+resolution,max_x,resolution*2)
        self.yrange = np.arange(min_y+resolution,max_y,resolution*2)

        self.xv, self.yv = np.meshgrid(self.xrange,self.yrange)

    def Create_KDTree(self):
        '''
        Init the kdtree
        '''
        return spatial.cKDTree(self.data.values) #input is an array

    def Count_Points(self):
        '''
        Compute the grid cell

        ROADMAP:
            Step 1. preselect points using KD-Tree
            Step 2. refine preselection using if statements (to get blocks)

        '''

        #time keeping
        start = time.time()

        radius = self.resolution*1.5
        grid_counter = np.zeros((len(self.yrange),len(self.xrange)))
        for idx_y,y in enumerate(self.yrange):
            for idx_x,x in enumerate(self.xrange):
                #step 1
                subset = self.kdtree.query_ball_point([x,y],r=radius)
                
                #step 2
                if not subset:
                    grid_counter[idx_y,idx_x] = 0
                else:
                    count = 0
                    for point in subset:
                        rdx,rdy = self.data.iloc[point].values
                        
                        if x-self.resolution < rdx and x+self.resolution > rdx and y-self.resolution < rdy and y+self.resolution > rdy:
                            count += 1
                    
                    grid_counter[idx_y,idx_x] = count
                
                lbrd = rijksdriehoek.Rijksdriehoek(x-self.resolution,y-self.resolution)
                rbrd = rijksdriehoek.Rijksdriehoek(x+self.resolution,y-self.resolution)
                rtrd = rijksdriehoek.Rijksdriehoek(x+self.resolution,y+self.resolution)
                ltrd = rijksdriehoek.Rijksdriehoek(x-self.resolution,y+self.resolution)

                leftbot = lbrd.to_wgs()
                rightbot = rbrd.to_wgs()
                righttop = rtrd.to_wgs()
                lefttop = ltrd.to_wgs()

                # leftbot = rijksdriehoek.rd_to_wgs(x-self.resolution,y-self.resolution)
                # rightbot = rijksdriehoek.rd_to_wgs(x+self.resolution,y-self.resolution)
                # righttop = rijksdriehoek.rd_to_wgs(x+self.resolution,y+self.resolution)
                # lefttop = rijksdriehoek.rd_to_wgs(x-self.resolution,y+self.resolution)

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

def main():
    '''
    main function to create the point maps

    '''
    print('starting timer for t37')
    start = time.time()

    filenamet37 = 't37_points.csv'
    objt37 = pointCounter(filenamet37)

    dft37 = objt37.start()
    dft37.to_csv(f't37_map_res200.csv')

    print('Done with track t37')
    print(f'Time elapsed is {time.time()-start} seconds')

    print('Starting with track t88')
    start = time.time()

    filenamet88 = 't88_points.csv'
    objt88 = pointCounter(filenamet88)

    dft88 = objt88.start()
    dft88.to_csv(f't88_map_res200.csv')

    print('Done with track t88')
    print(f'Time elapsed is {time.time()-start} seconds')

def _test():
    '''
    Testing of the class pointCounter
    '''
    obj = pointCounter()
    df = obj.start()
    df.to_csv('test3.csv')
    # print(obj.grid_counter)

    from mpl_toolkits.mplot3d import Axes3D  
    # Axes3D import has side effects, it enables using projection='3d' in add_subplot
    # import matplotlib.pyplot as plt
    # import random

    # def fun(x, y):
        # return x**2 + y

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    # x = y = np.arange(-3.0, 3.0, 0.05)
    # X, Y = np.meshgrid(x, y)
    # zs = np.array(fun(np.ravel(X), np.ravel(Y)))
    # Z = zs.reshape(X.shape)

    ax.plot_surface(obj.xv, obj.yv, obj.grid_counter,cmap='jet')

    ax.set_xlabel('X Label')
    ax.set_ylabel('Y Label')
    ax.set_zlabel('Z Label')

    # plt.show()

def _test2():
    '''
    testing conversion error
    '''
    FILENAME = 't37_points_d500.csv'
    obj = pointCounter(FILENAME)
    df = obj.start()

    # df = obj.Create_Points()
    df.to_csv('t37_point_subset_resolution.csv')

if __name__ == '__main__':
    main()
    #_test()
    # _test2()
