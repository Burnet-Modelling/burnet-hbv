# Generic utilities script
import pandas as pd
import openpyxl
import os
import socket
from platform import system

# Sharepoint folder, used for storing data and saving outputs
def _get_sharepoint_folder():
    """ Adapted from HCV code, function written by @kelly.maynard"""

    user = socket.gethostname()
    platform = system()

    if user in ['ChrisS-OCT25']:
        folder = r'C:\Users\chris.seaman\Burnet Institute\WG-Modelling-Hepatitis B - Documents\Applications\HepAus Submission Modelling'
    else:
        raise Exception(f'Error: unknown user "{user}", please add user information for future convenience!')

    return os.path.join(os.path.abspath(folder), '')

# GitHub repo, used for accessing frameworks, saving calibrations and databooks.
def _get_github_folder():
    """ Adapted from HCV code, function written by @kelly.maynard"""

    user = socket.gethostname()
    platform = system()

    if user in ['ChrisS-OCT25']:
        folder = r'C:\Users\chris.seaman\Desktop\GitRepos\burnet-hbv'
    else:
        raise Exception(f'Error: unknown user "{user}", please add user information for future convenience!')

    return os.path.join(os.path.abspath(folder), '')


# Dynamically read and import data from tables within an Excel workbook
def read_table(file_name: str, table_name: str) -> pd.DataFrame:

    """ Allows data to be read from tables in an excel workbook without needing to
    manually specify cell references.

    Taken from: https://stackoverflow.com/questions/54241345/pandas-read-a-table-from-excel"""

    wb = openpyxl.load_workbook(file_name, read_only= False, data_only = True) # openpyxl does not have table info if read_only is True; data_only means any functions will pull the last saved value instead of the formula
    for sheetname in wb.sheetnames: # pulls as strings
        sheet = wb[sheetname] # get the sheet object instead of string
        if table_name in sheet.tables: # tables are stored within sheets, not within the workbook, although table names are unique in a workbook
            tbl = sheet.tables[table_name] # get table object instead of string
            tbl_range = tbl.ref #something like 'C4:F9'
            break # we've got our table, bail from for-loop
    data = sheet[tbl_range] # returns a tuple that contains rows, where each row is a tuple containing cells
    content = [[cell.value for cell in row] for row in data] # loop through those row/cell tuples
    header = content[0] # first row is column headers
    rest = content[1:] # every row that isn't the first is data
    df = pd.DataFrame(rest, columns = header)
    wb.close()
    return df