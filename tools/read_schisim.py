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

# ----------------- Read gr3 file
def read_gr3(file_path: str) -> Tuple[pd.DataFrame, pd.DataFrame, List[List[str]]]:
    """
    Reads a GR3 file and extracts node and element data into pandas DataFrames, 
    along with any extra lines that don't fit standard node or element formats.

    Parameters:
    -----------
    file_path : str
        Path to the GR3 file.

    Returns:
    --------
    tuple:
        - nodes_df (pd.DataFrame): DataFrame containing nodes with columns ['id', 'long', 'lat', 'swe'].
        - elements_df (pd.DataFrame): DataFrame containing elements with columns 
          ['pg_id', 'num_nodes', 'node_id_1', 'node_id_2', 'node_id_3'].
        - extra_lines (List[List[str]]): List of lines that do not conform to standard node or element formats.

    Example:
    --------
    >>> nodes_df, elements_df, extra_lines = read_gr3("example.gr3")
    >>> print(nodes_df)
       id  long  lat  wse
    0   1   0.0  0.0  1.0
    1   2   1.0  0.0  2.0
    2   3   0.0  1.0  3.0
    >>> print(elements_df)
       pg_id  num_nodes  node_id_1  node_id_2  node_id_3
    0      1          3          1          2          3
    """

    nodes = []
    elements = []

    with open(file_path, 'r') as file:
        # Read the header line
        header = file.readline().strip()
        print("File Header:", header)

        node_lines = []
        element_lines = []
        extra_lines = []
        
        for line in file:
            line = line.strip().split()
            if len(line) == 4:  
                node_lines.append(line)
            elif len(line) == 5:  
                element_lines.append(line)
            elif len(line) > 5: 
                extra_lines.append(line)

        node_count = len(node_lines)
        element_count = len(element_lines)
        extra_count = len(extra_lines)
        print(f"Detected Nodes: {node_count}, Elements: {element_count}, Extra: {extra_count}")
        print(extra_lines)

        # Process node data
        for line in node_lines:
            node_id = int(line[0])
            x, y, depth = map(float, line[1:4])
            nodes.append({"id": node_id, "long": x, "lat": y, "wse": depth})

        # Process element data
        for line in element_lines:
            element_id = int(line[0])
            num_node, n1, n2, n3 = list(map(int, line[1:]))  # Nodes that form the element
            elements.append({
                "pg_id": element_id, 
                "num_nodes": num_node, 
                "node_id_1": n1, 
                "node_id_2": n2, 
                "node_id_3": n3
            })

    # Convert lists to DataFrames
    nodes_df = pd.DataFrame(nodes)
    elements_df = pd.DataFrame(elements)

    return nodes_df, elements_df, extra_lines

   