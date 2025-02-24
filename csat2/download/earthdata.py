import csat2
import sys
import os
import os.path
import logging
import json
import requests
from tqdm import tqdm

# For some reason, possibly relating to the Imperial network, IPv6 doesnt work well
# This might become an issue in the future, but it works for now.
requests.packages.urllib3.util.connection.HAS_IPV6 = False

log = logging.getLogger(__name__)

USERAGENT = "csat2/download_v{}".format(csat2.__version__).replace("\r", "")
TOKENFILE = os.environ["HOME"] + "/.csat2/laadsdaacrc"
AUTHFILE = os.environ["HOME"] + "/.csat2/earthdata_auth.json"
with open(AUTHFILE) as f:
    earthdata_auth = json.load(f)


# Based on code from
# https://urs.earthdata.nasa.gov/documentation/for_users/data_access/python
class SessionWithHeaderRedirection(requests.Session):
    AUTH_HOST = "urs.earthdata.nasa.gov"

    def __init__(self, username, password):
        super().__init__()
        self.auth = (username, password)

    # Overrides from the library to keep headers when redirected to or from
    # the NASA auth host.
    def rebuild_auth(self, prepared_request, response):
        headers = prepared_request.headers
        url = prepared_request.url

        if "Authorization" in headers:
            original_parsed = requests.utils.urlparse(response.request.url)
            redirect_parsed = requests.utils.urlparse(url)

            if (
                (original_parsed.hostname != redirect_parsed.hostname)
                and redirect_parsed.hostname != self.AUTH_HOST
                and original_parsed.hostname != self.AUTH_HOST
            ):
                del headers["Authorization"]
        return


def geturl(url, token=None, out=None, quiet=False):
    # session = requests.Session()
    session = SessionWithHeaderRedirection(
        earthdata_auth["username"], earthdata_auth["password"]
    )
    session.headers["user-agent"] = USERAGENT
    # if token is not None:
    #     session.headers["Authorization"] = f"Bearer {token}"
    log.debug("Setup session")
    response = session.get(url, stream=True)

    if response.status_code != 200:
        raise ValueError(f"Download failed status:{response.status_code}")
    total = int(response.headers.get("content-length", 0))
    log.debug("Got response")
    if out:
        with tqdm(
            desc=os.path.basename(url),
            total=total,
            unit="iB",
            unit_scale=True,
            unit_divisor=1024,
            disable=quiet,
        ) as bar:
            blocksize = max(4096, total // 100)
            for data in response.iter_content(chunk_size=blocksize):
                size = out.write(data)
                bar.update(size)
    else:
        return response.content


def get_token():
    try:
        with open(TOKENFILE) as f:
            return f.read().strip()
    except FileNotFoundError:
        raise IOError("Place your LAADS API key in the file {}".format(TOKENFILE))
