import logging
from .locator import configloader
import os

version = '0.1'

# Set logging level
logging.basicConfig(level=logging.ERROR)

# user_location - the default location for the user override file
user_location = os.path.expanduser('~/.csat2/config.cfg')

locator = configloader(user_location, default_fallback=True)
