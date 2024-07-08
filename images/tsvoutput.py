import csv
import argparse

def write_tsv(file_name, column_names, row_data):
    """
    Writes data to a TSV file.
    
    :param file_name: Name of the TSV file to write to.
    :param column_names: List of column names for the TSV file.
    :param row_data: List of row data to write to the TSV file.
    """
    data = [column_names, row_data]
    
    with open(file_name, 'w', newline='') as file:
        writer = csv.writer(file, delimiter='\t')
        writer.writerows(data)
    
    print(f"Data written to {file_name} successfully.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Create a TSV file with given data.')
    parser.add_argument('file_name', type=str, help='Name of the TSV file to create.')
    parser.add_argument('arg1', type=str, help='First argument (file path or sequencing ID).')
    parser.add_argument('arg2', type=str, help='Second argument (file path or sequencing ID).')
    
    args = parser.parse_args()
    
    column_names = ['individual_id', 'vcf']
    row_data = [args.arg1, args.arg2]
    
    write_tsv(args.file_name, column_names, row_data)
