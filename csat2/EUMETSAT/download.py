# %% Import EUMDAC and dependent libraries to begin
import eumdac
import datetime
import shutil
import requests
# %%

# Import EUMDAC and dependent libraries to begin

token = eumdac.AccessToken(("KEY", "TOKEN"))

datastore = eumdac.DataStore(token)

high_res = 'EO:EUM:DAT:0662'
normal_res = 'EO:EUM:DAT:0661'

selected_product = datastore.get_product(
    product_id='W_XX-EUMETSAT-Darmstadt,IMG+SAT,MTI1+FCI-1C-RRAD-FDHSI-FD--x-x---x_C_EUMT_20250115090301_IDPFI_OPE_20250115090007_20250115090924_N__O_0055_0000',
    collection_id='EO:EUM:DAT:0662')

try:
    with selected_product.open() as fsrc, \
            open(fsrc.name, mode='wb') as fdst:
        shutil.copyfileobj(fsrc, fdst)
    print(f'Download of product {selected_product} finished.')
except eumdac.product.ProductError as error:
    print(f"Error related to the product '{selected_product}' while trying to download it: '{error.msg}'")
except requests.exceptions.ConnectionError as error:
    print(f"Error related to the connection: '{error.msg}'")
except requests.exceptions.RequestException as error:
    print(f"Unexpected error: {error}")
# %%
