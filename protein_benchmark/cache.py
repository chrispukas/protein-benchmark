import os
import pickle

from importlib.resources import files



class Cache():
    def __init__(self, cache_dir=None, pickle_name=None):
        if cache_dir is None:
            self.cache_dir = os.path.join(os.path.dirname(files('protein_benchmark')), 'cache')
        else:
            self.cache_dir = cache_dir
        
        self.pickle_name = pickle_name

        print("Initialized cache directory:", self.cache_dir)

        if not os.path.exists(self.cache_dir):
            os.makedirs(self.cache_dir, exist_ok=True)

    def load_pickle(self, pickle_name=None) -> object:
        """Load a pickle file."""
        file_path = os.path.join(self.cache_dir, f"{self.pickle_name if pickle_name is None else pickle_name}.pkl")
        with open(file_path, 'rb') as file:
            return pickle.load(file, encoding='latin1')
    
    def save_pickle(self, data, pickle_name=None):
        """Save data to a pickle file."""
        file_path = os.path.join(self.cache_dir, f"{self.pickle_name if pickle_name is None else pickle_name}.pkl")
        with open(file_path, 'wb') as file:
            pickle.dump(data, file)
            
    def is_pickle(self, pickle_name=None) -> bool:
        """Check if a pickle file exists."""
        file_path = os.path.join(self.cache_dir, f"{self.pickle_name if pickle_name is None else self.pickle_name}.pkl")
        res = os.path.exists(file_path)
        print("Checking if pickle exists at:", file_path, "->", res)
        return res
    
def c_template(cache_name, compute_fn, *args, **kwargs):
    cache = Cache(pickle_name=cache_name)
    if cache_name and cache.is_pickle():
        return cache.load_pickle()
    else:
        result = compute_fn(*args, **kwargs)
        cache.save_pickle(result)
        return result