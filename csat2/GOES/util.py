_resolutions = [1, 0.5, 1, 2, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2]


def get_resolution(channel):
    """Returns the band resolution in km at nadir"""
    return _resolutions[channel - 1]
