# tools/read_schisim.py
import pandas as pd
from typing import List, Tuple

# ----------------- Read 2dm file
def read_2dm(file_path: str) -> Tuple[List[List[float]], List[List[int]], pd.DataFrame, pd.DataFrame]:
    """
    Reads a 2DM file and extracts node and element data into structured lists and pandas DataFrames.

    Parameters:
    -----------
    file_path : str
        Path to the 2DM file.

    Returns:
    --------
    tuple:
        - node_data (List[List[float]]): A list of nodes where each node is represented as [node_id, x, y, z].
        - element_data (List[List[int]]): A list of elements where each element is represented as [element_id, node1, node2, node3].
        - nodes_df (pd.DataFrame): DataFrame containing nodes with columns ['node_id', 'x', 'y', 'z'].
        - elements_df (pd.DataFrame): DataFrame containing elements with columns ['element_id', 'node1', 'node2', 'node3'].

    Example:
    --------
    >>> node_data, element_data, nodes_df, elements_df = read_2dm("example.2dm")
    >>> print(nodes_df)
       node_id    x    y    z
    0        1  0.0  0.0  0.0
    1        2  1.0  0.0  0.0
    2        3  0.0  1.0  0.0
    >>> print(elements_df)
       element_id  node1  node2  node3
    0           1      1      2      3
    """
    node_data = []
    element_data = []

    with open(file_path, 'r') as file:
        for line in file:
            tokens = line.split()
            if not tokens:
                continue
            # Parse node entries (ND) and store node ID with x, y, z coordinates
            if tokens[0] == 'ND':
                node_id = int(tokens[1])
                x, y, z = map(float, tokens[2:5])
                node_data.append([node_id, x, y, z])
            # Parse triangular element entries (E3T) with 4 columns (ID + 3 node indices)
            elif tokens[0] == 'E3T':
                element_id = int(tokens[1])
                node_indices = list(map(int, tokens[2:5]))  # 3 node indices
                element_data.append([element_id] + node_indices)

    # Convert to pandas DataFrames
    nodes_df = pd.DataFrame(node_data, columns=['node_id', 'x', 'y', 'z'])
    elements_df = pd.DataFrame(element_data, columns=['element_id', 'node1', 'node2', 'node3'])

    return node_data, element_data, nodes_df, elements_df

