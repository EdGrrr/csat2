from csat2 import locator
from .download import download_file_locations
import os


def create_metadata_year(year):
    output_folder = locator.get_folder("CALIPSO", "GEOMETA", year=year)

    try:
        os.makedirs(output_folder)
    except FileExistsError:
        pass

    outputfile = f"{output_folder}/calipso_{year}.txt"

    product = "LID_L1"
    col = "V4-11"

    filenames = []
    for mon in range(1, 13):
        try:
            nfilenames = download_file_locations(product, year, mon, col=col)
            filenames.extend(nfilenames)
        except:
            pass

    with open(outputfile, "w") as f:
        f.write("\n".join(filenames))


def create_metadata():
    for year in range(2006, 2022):
        create_metadata_year(year)
        print(f"{year} complete")
