
import numpy as np

class Data:
    def __init__(self, x_arr, y_arr):
        self.x = np.array(x_arr)
        self.y = np.array(y_arr)
        self.log_x = np.log10(abs(x_arr))
        self.log_y = np.log10(abs(y_arr))        
        
                


    
