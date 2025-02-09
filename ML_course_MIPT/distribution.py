import numpy as np

class LaplaceDistribution:    
    @staticmethod
    def mean_abs_deviation_from_median(x: np.ndarray):
        '''
        Args:
        - x: A numpy array of shape (n_objects, n_features) containing the data
          consisting of num_train samples each of dimension D.
        '''
        ####
        # Do not change the class outside of this block
        mediana = np.median(x, axis=0)
        return np.mean(abs(np.subtract(x, mediana)), axis=0)
        # Your code here
        ####

    def __init__(self, features):
        '''
        Args:
            feature: A numpy array of shape (n_objects, n_features). Every column represents all available values for the selected feature.
        '''
        ####
        # Do not change the class outside of this block
        self.loc = np.median(features, axis=0) # YOUR CODE HERE мю
        self.scale = self.mean_abs_deviation_from_median(features)# YOUR CODE HERE b
        ####


    def logpdf(self, values):
        '''
        Returns logarithm of probability density at every input value.
        Args:
            values: A numpy array of shape (n_objects, n_features). Every column represents all available values for the selected feature.
        '''
        ####
        # Do not change the class outside of this block
        print(np.emath.log(1/(2 * self.scale) * np.exp(-abs(np.subtract(values, self.loc))/self.scale)))
        return np.emath.log(1/(2 * self.scale) * np.exp(-abs(np.subtract(values, self.loc)/self.scale)))
        ####
        
    
    def pdf(self, values):
        '''
        Returns probability density at every input value.
        Args:
            values: A numpy array of shape (n_objects, n_features). Every column represents all available values for the selected feature.
        '''
        return np.exp(self.logpdf(values))
