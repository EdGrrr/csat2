from ._version import __version__
from .locator import configloader
import os

# user_location - the default location for the user override file
user_location = os.path.expanduser('~/.csat2/config.cfg')

locator = configloader(user_location, default_fallback=True)
