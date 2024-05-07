import os
import shutil
import logging
import time
import csv
from functools import wraps


# Setup logger
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)
formatter = logging.Formatter('%(levelname)s:%(asctime)s:%(name)s:%(message)s')

# Define handlers
fh = logging.FileHandler('{}.log'.format(__name__))
fh.setFormatter(formatter)
sh = logging.StreamHandler()
sh.setLevel(logging.WARNING)

# Add handlers to logger
logger.addHandler(fh)
logger.addHandler(sh)


def timer(orig_func):
    @wraps(orig_func)
    def wrapper(*args, **kwargs):
        t1 = time.time()
        result =  orig_func(*args, **kwargs)
        t2 = time.time() - t1
        logger.info('{} ran in {} sec'.format(orig_func.__name__, round(t2, 2)))

        # Uncomment below to record the timer results
        # with open('timer.csv', 'a') as file:
        #     writer = csv.writer(file, delimiter=',')
        #     writer.writerows([(str(orig_func.__name__), t2)])

        return result
    return wrapper


def check_dirpath(dirpath) -> str:
    """Check if input dir paths ends with '/', if not concat '/'."""
    return dirpath + '/' if not dirpath.endswith('/') else dirpath


def create_dir(tail_dir, head_dir='', replace=False, verbose=False) -> str:
    """Create a new dir from input tail and head. If already exists, remove existing dir first. Note: only works for UNIX systems.

    Args:
        tail_dir (str): dirpath to be created.
        head_dir (str, optional): a pre-existing dirpath. Defaults to ''.
        replace (bool, optional): if True, replace if already exists. Defaults to True.
        verbose (bool, optional): Defaults to False.

    Returns:
        str: created dirpath.
    """

    # Ensure dirpaths end with a '/'
    ends_with = lambda dirpath: dirpath + '/' if not dirpath.endswith('/') and dirpath != '' else dirpath

    # Concat tail and head dirpaths
    new_dir = ends_with(tail_dir) + head_dir

    # Remove new dir if it already exists, else make new dir
    if os.path.exists(new_dir) and os.path.isdir(new_dir):
        if replace:
            shutil.rmtree(new_dir)
            os.makedirs(new_dir)
            logger.info('Warning! Replaced old {}'.format(new_dir))
        else:
            raise FileExistsError('{} directory already exists. To replace old directory, set replace=True.'.format(new_dir))
    else:
        os.makedirs(new_dir)
        logger.info('Created {}'.format(new_dir))

    return new_dir
