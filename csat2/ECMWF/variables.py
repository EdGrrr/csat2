# This translates between local names and the cdsapi names
# Partly to ensure independence if they decide to change names
# in the future
# The list is ['local_name', 'cds_name']
CDS_NAME_TABLE = [
    ["U-wind-component", "u_component_of_wind"],
    ["V-wind-component", "v_component_of_wind"],
    ["10m-U-wind-component", "10m_u_component_of_wind"],
    ["10m-V-wind-component", "10m_v_component_of_wind"],
    ["10m-U-wind-component", "10_metre_u_wind_component"],
    ["10m-V-wind-component", "10_metre_v_wind_component"],
    ["SST", "sea_surface_temperature"],
    ["Mean_sea_level_pressure", "mean_sea_level_pressure"],
]

LOCAL_NAMES, CDS_NAMES = list(zip(*CDS_NAME_TABLE))


def _convert_vname_to_cds(vname):
    try:
        return CDS_NAMES[LOCAL_NAMES.index(vname)]
    except:
        return vname.lower()


def _convert_cds_to_vname(cds_name):
    try:
        return LOCAL_NAMES[CDS_NAMES.index(cds_name)]
    except:
        return cds_name.capitalize()


def _get_levelstr(level):
    if level == "surf":
        levelstr = "surf"
    else:  # level is some kind of pressure level
        if isinstance(level, str):
            level = level.replace("hPa", "")
        if isinstance(level, int):
            level = "{.0f}".format(level)
        levelstr = "{}hPa".format(level)
    return levelstr