import io
import os
import zipfile
from knool.satellite.data import packet_decoding
import pandas as pd
import numpy as np
from ...fortlib import fdbaqdecoder
from ..algorithm.sentinel1decoder.sentinel1decoder import _headers as hdrs
from ..algorithm.sentinel1decoder.sentinel1decoder import constants as cnst
from ..algorithm.sentinel1decoder.sentinel1decoder.utilities import read_subcommed_data
from ..algorithm.sentinel1decoder.sentinel1decoder import utilities


class Sentinel_L0IW():

    def __init__(self, zip_file):
        self.zip_file = zip_file
        self.dat_files = []
        with zipfile.ZipFile(zip_file, 'r') as zip_data:
            infos = zip_data.infolist()
            for i, info in enumerate(infos):
                _, ext = os.path.splitext(info.filename)
                if ext == ".dat":
                    self.dat_files.append(info.filename)
                    kind = os.path.splitext(info.filename)[0][-5:]
                    if (kind != "index") and (kind != "annot"):
                        self.bytesIO = io.BytesIO(zip_data.read(info.filename))

        # Split metadata into blocks of consecutive packets w/ const swath number
        self._packet_metadata = self._index_df_on_bursts(
            self.decode_metadata())
        self._ephemeris = None
        self._burst_data_dict = dict.fromkeys(
            self._packet_metadata.index.unique(
                level=cnst.BURST_NUM_FIELD_NAME),
            None
        )

    @property
    def packet_metadata(self) -> pd.DataFrame:
        return self._packet_metadata

    @property
    def ephemeris(self) -> pd.DataFrame:
        if self._ephemeris is None:
            self._ephemeris = read_subcommed_data(self.packet_metadata)
        return self._ephemeris

    def get_burst_metadata(self, burst: int) -> pd.DataFrame:
        return self.packet_metadata.loc[burst]

    def get_burst_data(self, burst: int, huffman_tree) -> np.array:
        self._burst_data_dict[burst] = self.decode_packets(
            self.get_burst_metadata(burst), huffman_tree)
        return self._burst_data_dict[burst]

    def get_huffman_tree(self):
        _TREE_BRC_ZERO = (0, (1, (2, 3)))
        _TREE_BRC_ONE = (0, (1, (2, (3, 4))))
        _TREE_BRC_TWO = (0, (1, (2, (3, (4, (5, 6))))))
        _TREE_BRC_THREE = ((0, 1), (2, (3, (4, (5, (6, (7, (8, 9))))))))
        _TREE_BRC_FOUR = (
            (0, (1, 2)), ((3, 4), ((5, 6), (7, (8, (9, ((10, 11), ((12, 13), (14, 15)))))))))

        huffman_tree = []
        huffman_tree.append(packet_decoding.build_huffman_tree(_TREE_BRC_ZERO))
        huffman_tree.append(packet_decoding.build_huffman_tree(_TREE_BRC_ONE))
        huffman_tree.append(packet_decoding.build_huffman_tree(_TREE_BRC_TWO))
        huffman_tree.append(
            packet_decoding.build_huffman_tree(_TREE_BRC_THREE))
        huffman_tree.append(packet_decoding.build_huffman_tree(_TREE_BRC_FOUR))
        huffman_tree = packet_decoding.convert_huffman_tree_fortran(
            huffman_tree).astype(np.int32)
        return huffman_tree

    def decode_metadata(self) -> pd.DataFrame:
        output_row_list = []
        self.bytesIO.seek(0)
        while True:
            try:
                pos = self.bytesIO.tell()
                output_dictionary_row, _ = self._read_single_packet(
                    self.bytesIO)
            except:
                break
            output_dictionary_row["Byte_Position"] = pos
            output_row_list.append(output_dictionary_row)

        output_dataframe = pd.DataFrame(output_row_list)
        self.bytesIO.seek(0)
        return output_dataframe

    def decode_packets(self, input_header: pd.DataFrame, huffman_tree) -> np.array:
        self.bytesIO.seek(input_header["Byte_Position"].min())
        swath_numbers = input_header[cnst.SWATH_NUM_FIELD_NAME].unique()
        num_quads = input_header[cnst.NUM_QUADS_FIELD_NAME].unique()

        packet_counter = 0
        packets_to_process = len(input_header)
        nq = input_header[cnst.NUM_QUADS_FIELD_NAME].unique()[0]
        print(input_header["Byte_Position"].min(),
              nq, num_quads, swath_numbers)

        max_pdl = input_header["Packet Data Length"].max()
        num_baq_blocks = np.ceil(nq/128)
        num_bytes_list = []
        output_data = np.zeros([packets_to_process, nq * 2], dtype=(complex))
        byte_array = np.zeros([max_pdl, packets_to_process], dtype=np.uint8)
        while packet_counter < packets_to_process:
            this_header, packet_data_bytes = self._read_single_packet(
                self.bytesIO)
            test = np.frombuffer(packet_data_bytes, dtype=np.uint8)
            num_bytes_list.append(test.shape[0])
            byte_array[0:num_bytes_list[packet_counter], packet_counter] = test
            packet_counter += 1
        param_dir = "/home/knakata/program/mycode/kn-general_analysis_tools/fortlib/share/sentinel1_decorder/"
        print(byte_array.shape, num_baq_blocks, nq, np.array(
            num_bytes_list).shape, huffman_tree.shape)
        output_data = fdbaqdecoder.allproc(byte_array, num_baq_blocks, nq, np.array(
            num_bytes_list), huffman_tree, param_dir)
        self.bytesIO.seek(0)
        return output_data.T

    def _generate_burst_cache_filename(self, burst: int) -> str:
        return os.path.splitext(self.filename)[0] + "_b" + str(burst) + ".npy"

    def _index_df_on_bursts(self, packet_metadata: pd.DataFrame) -> pd.DataFrame:
        packet_metadata = packet_metadata.groupby(
            packet_metadata[[cnst.SWATH_NUM_FIELD_NAME,
                             cnst.NUM_QUADS_FIELD_NAME]]
            .diff()
            .ne(0)
            .any(axis=1)
            .cumsum(), group_keys=True
        ).apply(lambda x: x)

        packet_metadata.index.names = [
            cnst.BURST_NUM_FIELD_NAME,
            cnst.PACKET_NUM_FIELD_NAME,
        ]

        return packet_metadata

    def _read_single_packet(self, bytesIO):
        data_buffer = bytesIO.read(6)
        output_dictionary_row = hdrs.decode_primary_header(data_buffer)
        pkt_data_len = output_dictionary_row[cnst.PACKET_DATA_LEN_FIELD_NAME]
        packet_data_buffer = bytesIO.read(pkt_data_len)
        secondary_hdr = hdrs.decode_secondary_header(packet_data_buffer[:62])
        output_dictionary_row.update(secondary_hdr)
        output_bytes = packet_data_buffer[62:]
        return output_dictionary_row, output_bytes
