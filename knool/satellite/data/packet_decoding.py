import numpy as np


def build_huffman_tree(huffman_list):
    """
    ネストされたリストを Huffman 木の辞書形式に変換する。

    :param huffman_list: Huffman 木のネストされたタプル
    :return: Huffman 木の辞書
    """
    huffman_tree = {}
    node_counter = [0]  # ノードIDのカウンター（リストでミュータブルにする）

    def add_node(subtree):
        """再帰的にツリーを追加する"""
        node_id = node_counter[0]  # 現在のノードID
        node_counter[0] += 1  # ノードIDを更新

        if isinstance(subtree, tuple):  # ブランチノード（リスト・タプル）
            left_child = add_node(subtree[0])  # 左の子
            right_child = add_node(subtree[1])  # 右の子
            huffman_tree[node_id] = (left_child, right_child)  # 子ノードを登録
        else:  # リーフノード（文字列・数値）
            huffman_tree[node_id] = subtree

        return node_id  # 親ノードにノードIDを返す

    add_node(huffman_list)  # ルートノードを構築
    return huffman_tree


def convert_huffman_tree_fortran(huffman_tree_list):
    num_type = len(huffman_tree_list)
    num_node = 0
    for huffman_tree in huffman_tree_list:
        num_node = np.max([num_node, len(huffman_tree)])

    huffman_tree_fortran = np.zeros([num_type, num_node, 2])
    huffman_tree_fortran[:, :, :] = -999
    for i, huffman_tree in enumerate(huffman_tree_list):
        for node in range(len(huffman_tree)):
            if isinstance(huffman_tree[node], int):
                huffman_tree_fortran[i, node, :] = -huffman_tree[node]
            else:
                huffman_tree_fortran[i, node, 0] = huffman_tree[node][0]
                huffman_tree_fortran[i, node, 1] = huffman_tree[node][1]
    return huffman_tree_fortran
