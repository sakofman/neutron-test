import csv

"""
File processing module
"""

# TODO: add functionality for converting a string (if we have to specify an exponential or smth, like for the beam parameters) to a float if encountered 
def params_to_dict(filepath, delimiter=None):
    """
    Returns a dict containing object parameters specified in the input file.
    The text file should have the first row as keys and the second row as values.
    If no delimiter is provided, the function will split on any whitespace.

    Arguments
    ---------
    filepath : string
        Path to file.
    delimiter : string (optional)
        Delimiter to use for splitting the keys and values; default is a space.
        
    Returns
    -------
    dict
        Dictionary with keys as the first line of the file read and corresponding 
        values as the columns below.
    """
    params_dict = {}

    with open(filepath, mode='r') as file:
        lines = file.readlines()
        # Remove leading/trailing whitespace and ignore blank lines
        non_blank_lines = [line.strip() for line in lines if line.strip()]

    if len(non_blank_lines) != 2:
        print("Warning: The file should contain exactly two lines of data. Using the first line as keys and the second line as values.")
    
    # Split lines using provided delimiter
    if delimiter:
        keys = non_blank_lines[0].split(delimiter)
        values = non_blank_lines[1].split(delimiter)
    else:
        keys = non_blank_lines[0].split()
        values = non_blank_lines[1].split()

    # Combine read data in a dict
    for key, value in zip(keys, values):
        params_dict[key] = float(value)

    return params_dict
