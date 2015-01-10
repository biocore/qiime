#!/usr/bin/env python
"""tests for flowgram clustering"""

__author__ = "Jens Reeder"
__copyright__ = "Copyright 2011, The QIIME Project"
# remember to add yourself if you make changes
__credits__ = ["Jens Reeder", "Rob Knight", "Nigel Cook"]
__license__ = "GPL"
__version__ = "1.9.0-rc2"
__maintainer__ = "Jens Reeder"
__email__ = "jens.reeder@gmail.com"

import signal
import os
from os import rmdir
from time import sleep
from StringIO import StringIO
from tempfile import mkdtemp

from unittest import TestCase, main
from numpy.testing import assert_almost_equal
from bfillings.denoiser import Flowgram, parse_sff, FlowgramCollection
from burrito.util import ApplicationNotFoundError

from qiime.denoiser.flowgram_clustering import *
from qiime.denoiser.utils import FlowgramContainerArray
from qiime.denoiser.cluster_utils import setup_workers, setup_server, stop_workers

# timeout handling taken from test_workflow.py


class TimeExceededError(Exception):
    pass

allowed_seconds_per_test = 60


def timeout(signum, frame):
    raise TimeExceededError("Test failed to run in allowed time (%d seconds)."
                            % allowed_seconds_per_test)


class DenoiserTests(TestCase):

    def setUp(self):
        self.rec = data
        self.flowgram = flowgram
        self.flowgrams, self.head = parse_sff(self.rec.split("\n"))
        self.socket = None
        self.tmp_dir = None

        signal.signal(signal.SIGALRM, timeout)
        # set the 'alarm' to go off in allowed_seconds seconds
        signal.alarm(allowed_seconds_per_test)

    def tearDown(self):

        # turn off the alarm
        signal.alarm(0)

        if self.socket:
            self.socket.close()
        if self.tmp_dir:
            try:
                rmdir(self.tmp_dir)
            except OSError:
                # give client some more time to shut down and clean up
                sleep(3)
                rmdir(self.tmp_dir)

    def test_compute_workload(self):
        """allocation of workload based on observed worker performance"""

        spread = [1.23, 0.5, 1.27]
        num_cores = 3
        num_flows = 11
        result = compute_workload(num_cores, num_flows, spread)
        self.assertEqual(result, [4, 2, 5])

    def test_adjust_processing_time(self):
        """evalulate worker performance"""

        workload = [6, 10, 16]
        num_cores = 3
        timing = [9.0, 11.0, 65.0]
        epoch = 1.0
        result = adjust_processing_time(num_cores, workload, timing, epoch)
        self.assertEqual(result, [1.125, 1.5, 0.375])

    def test_get_flowgram_distances(self):
        """get_flowgram_distances compute the correct alignment score."""

        scores, names, fc = get_flowgram_distances(
            "1", self.flowgram, self.flowgrams,
            FlowgramContainerArray(),
            {"FZTHQMS01CIW5N": ""}, "/tmp/")
        self.assertEqual(names, ["FZTHQMS01CIW5N"])
        assert_almost_equal(scores, [[4.95274923, 0.7815385]], decimal=4)

    def test_get_flowgram_distances_on_cluster(self):
        """get_flowgram_distances_on_cluster computes the correct alignment score."""

        self.tmp_dir = mkdtemp(dir="./", suffix="/")
        # setup server and workers
        self.socket = setup_server()
        workers, client_sockets = setup_workers(1, self.tmp_dir, self.socket,
                                                verbose=False)
        client_sockets = [a for a, b in client_sockets]
        scores, names, fc = get_flowgram_distances_on_cluster(
            "1", self.flowgram,
            self.flowgrams,
            FlowgramContainerArray(
            ),
            {"FZTHQMS01CIW5N": ""},
            1, 1, [1],
            client_sockets)
        stop_workers(client_sockets)

        self.assertEqual(names, ["FZTHQMS01CIW5N"])
        assert_almost_equal(scores, [[4.95274923, 0.7815385]], decimal=4)

    def test_log_remaining_rounds(self):
        """We can calculate how far we have to go"""

        # all empty
        self.assertEqual(log_remaining_rounds(dict(), dict(), 0), 0)

        # with something in it
        ids = dict({'1': 1, '2': 1, '3': 1, '4': 1})
        mapping = dict({'1': [5, 6], '2': [], '3': [7], '4': [8, 9, 10]})

        self.assertEqual(log_remaining_rounds(ids, mapping, 0), 4)
        self.assertEqual(log_remaining_rounds(ids, mapping, 1), 3)
        self.assertEqual(log_remaining_rounds(ids, mapping, 2), 2)
        self.assertEqual(log_remaining_rounds(ids, mapping, 5), 0)

flowgram = Flowgram("0.99	0.00	0.99	0.00	0.00	1.02	0.00	1.00	1.00	1.12\t\t" +
                    "0.01	0.01	1.89	0.01	0.95	0.95	0.97	0.00	0.02	0.98\t" +
                    "0.97	0.00	0.97	0.05	0.01	1.06	0.03	0.97	0.00	0.03\t" +
                    "0.97	0.02	0.00	1.09	0.02	0.01	0.96	0.00	0.00	1.01\t" +
                    "0.04	0.00	0.99	0.06	0.97	0.00	0.09	0.97	0.04	0.00\t" +
                    "1.94	0.09	1.02	0.00	2.86	1.02	1.00	1.11	0.10	1.97\t" +
                    "0.12	0.98	0.01	0.99	2.90	0.03	0.04	1.93	0.15	1.02\t" +
                    "1.95	1.00	1.02	0.00	0.12	1.00	0.97	0.00	1.00	0.06\t" +
                    "0.97	0.00	0.96	0.05	0.10	1.03	0.12	0.99	1.98	0.09\t" +
                    "1.99	0.08	0.13	2.10	0.14	0.05	1.00	0.10	0.00	1.00\t" +
                    "1.00	0.00	0.07	4.82	0.10	1.04	2.05	0.00	2.01	0.04\t" +
                    "1.96	0.08	0.93	0.00	0.93	0.03	0.99	0.02	1.01	0.06\t" +
                    "0.09	1.04	0.14	1.06	0.07	2.04	3.49	0.15	1.02	0.80\t" +
                    "0.23	0.07	1.07	0.17	1.91	0.07	0.18	1.00	0.32	0.07\t" +
                    "0.97	0.11	0.96	0.96	0.14	1.96	0.19	2.01	2.84	0.28\t" +
                    "0.08	2.03	1.32	0.06	0.05	1.10	0.17	0.88	0.09	0.95\t" +
                    "0.14	0.13	1.85	1.07	1.78	0.89	1.94	0.19	1.09	0.14\t" +
                    "1.09	0.13	0.13	0.86	1.85	0.07	0.09	1.97	1.20	0.08\t" +
                    "0.95	0.23	0.09	0.94	0.16	0.11	1.92	0.12	0.89	1.95\t" +
                    "0.21	0.12	0.97	0.14	0.16	1.86	0.12	1.89	1.00	1.07\t" +
                    "0.06	0.16	1.05	0.11	0.06	0.95	0.12	0.13	1.01	0.15\t" +
                    "3.79	0.14	0.15	0.98	0.40	0.11	1.00	0.19	1.01	1.09\t" +
                    "0.12	0.94	0.11	0.15	1.00	2.04	2.03	0.95	0.06	3.05\t" +
                    "0.22	0.08	1.82	0.21	1.02	0.09	2.88	1.88	0.15	0.07\t" +
                    "1.05	1.89	0.08	0.06	1.87	2.87	1.87	0.06	0.15	1.15\t" +
                    "0.25	0.08	0.96	0.12	0.06	0.95	0.09	0.13	1.05	1.95\t" +
                    "3.81	1.02	0.13	0.17	2.14	1.08	0.19	0.13	1.08	1.01\t" +
                    "1.99	0.11	0.18	1.06	0.17	0.04	0.98	0.08	1.01	2.86\t" +
                    "1.06	0.96	0.10	0.22	1.99	2.04	0.14	0.00	0.97	0.16\t" +
                    "0.95	0.07	2.75	0.02	0.98	0.12	2.94	0.00	0.99	1.03\t" +
                    "0.26	2.89	0.15	1.87	0.10	0.15	0.98	0.17	1.07	0.92\t" +
                    "0.00	0.09	1.08	0.16	3.78	1.01	0.07	0.87	0.22	0.98\t" +
                    "1.97	1.09	0.08	0.17	1.08	0.03	0.97	2.04	0.18	0.14\t" +
                    "1.03	0.03	0.00	1.16	0.12	1.81	2.06	0.18	0.17	2.06\t" +
                    "0.14	0.85	0.21	0.12	1.01	1.05	1.05	0.94	0.99	0.11\t" +
                    "0.15	1.08	2.00	1.02	0.99	0.13	1.07	0.13	0.98	0.16\t" +
                    "0.09	0.99	3.00	1.05	1.02	0.02	0.10	0.93	0.11	0.09\t" +
                    "0.81	0.97	0.13	0.05	2.04	1.93	1.12	0.04	0.93	0.93\t" +
                    "0.11	0.06	1.96	0.06	0.09	1.14	0.15	0.06	1.08	0.06\t" +
                    "0.94	0.11	0.00	0.88	1.11	0.10	2.08	1.05	0.15	0.09")

data = """Common Header:
  Magic Number:  0x2E736666
  Version:       0001
  Index Offset:  7773224
  Index Length:  93365
  # of Reads:    1
  Header Length: 440
  Key Length:    4
  # of Flows:    400
  Flowgram Code: 1
  Flow Chars:    TACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACG
  Key Sequence:  TCAG

>FZTHQMS01CIW5N

  Run Prefix:   R_2009_07_27_13_43_24_
  Region #:     1
  XY Location:  0918_1993

  Run Name:       R_2009_07_27_13_43_24_FLX04080350_Administrator_MIS_7_27_09
  Analysis Name:  D_2009_07_27_17_44_53_FLX04080350_fullProcessingAmplicons
  Full Path:      /data/R_2009_07_27_13_43_24_FLX04080350_Administrator_MIS_7_27_09/D_2009_07_27_17_44_53_FLX04080350_fullProcessingAmplicons/

  Read Header Len:  32
  Name Length:      14
  # of Bases:       271
  Clip Qual Left:   17
  Clip Qual Right:  176
  Clip Adap Left:   0
  Clip Adap Right:  0

Flowgram:	0.99	0.00	0.97	0.00	0.01	1.04	0.00	0.96	1.01	1.08	0.00	0.00	1.99	0.00	0.95	1.00	1.01	0.00	0.00	1.00	1.03	0.00	0.96	0.03	0.01	1.00	0.00	1.00	0.01	0.00	0.96	0.01	0.00	1.07	0.00	0.00	0.98	0.00	0.00	0.99	0.03	0.00	0.97	0.01	0.97	0.00	0.01	1.02	0.05	0.00	2.03	0.04	1.02	0.00	3.08	1.01	1.03	1.04	0.07	2.06	0.11	0.99	0.00	0.96	3.03	0.05	0.02	2.01	0.15	1.02	2.05	1.00	0.99	0.00	0.04	1.01	1.00	0.00	0.97	0.05	1.00	0.00	0.98	0.05	0.14	0.96	0.06	1.02	2.05	0.09	2.02	0.06	0.14	2.00	0.10	0.05	1.01	0.07	0.00	1.01	1.01	0.03	0.04	4.99	0.09	1.02	2.07	0.02	2.02	0.03	2.05	0.07	0.97	0.00	0.99	0.03	0.98	0.05	0.98	0.10	0.08	1.03	0.04	1.01	0.13	1.98	4.17	0.12	1.05	0.96	0.10	0.00	1.00	0.12	2.04	0.03	0.15	0.98	0.04	0.08	1.02	0.12	1.04	1.03	0.02	2.01	0.12	2.01	3.03	0.00	0.00	1.98	0.96	0.08	0.01	0.96	0.16	1.03	0.02	1.01	0.08	0.13	2.24	1.04	2.11	1.08	2.12	0.12	0.97	0.06	0.93	0.05	0.11	0.97	1.95	0.05	0.10	2.02	0.95	0.11	1.05	0.00	0.13	1.06	0.10	0.11	2.05	0.14	1.02	1.99	0.00	0.01	0.93	0.08	0.11	2.04	0.16	1.99	1.09	0.93	0.05	0.00	1.00	0.13	0.12	0.94	0.00	0.02	1.04	0.04	4.01	0.08	0.18	0.94	0.00	0.08	1.05	0.08	1.06	1.05	0.09	0.93	0.00	0.04	0.99	1.85	2.02	0.98	0.00	2.69	0.10	0.12	1.93	0.14	1.04	0.01	3.00	2.06	0.13	0.00	0.79	1.98	0.50	0.06	1.01	2.47	1.08	0.59	0.11	0.97	1.04	0.12	0.54	0.45	0.03	0.49	0.47	0.08	0.51	1.50	2.13	0.54	0.53	1.00	1.54	0.54	0.50	0.01	1.09	0.53	1.00	0.59	0.04	0.92	0.09	0.13	1.63	0.51	0.55	1.51	0.97	0.92	1.07	0.00	1.01	1.51	0.13	0.03	0.99	0.09	1.04	1.56	1.94	0.52	0.38	0.00	2.46	0.90	0.51	0.48	0.43	1.30	0.51	0.96	1.50	0.08	0.92	0.07	2.00	0.52	0.42	0.50	0.43	1.52	1.96	1.60	0.09	0.29	0.51	0.38	1.64	1.04	0.01	0.17	0.98	0.05	2.86	1.74	0.09	0.45	0.31	0.45	1.13	0.97	0.13	1.01	1.52	0.02	0.66	2.13	0.11	0.52	0.46	0.18	0.63	0.94	0.45	1.68	1.64	0.11	0.08	1.47	0.95	1.08	0.64	0.01	1.04	0.43	0.98	0.61	0.60	0.41	1.49	0.99	1.54	0.41	0.43	0.47	0.38	0.11	1.01	0.52	0.19	0.41	2.57	1.55	1.03	0.07	0.26	0.94	0.03	0.04	1.57	0.50	0.10	0.43	0.88	1.10	1.08	0.12	0.87	0.35	0.18	0.45	1.63	0.08	1.04	0.93	0.00	0.15
Flow Indexes:	1	3	6	8	9	10	13	13	15	16	17	20	21	23	26	28	31	34	37	40	43	45	48	51	51	53	55	55	55	56	57	58	60	60	62	64	65	65	65	68	68	70	71	71	72	73	76	77	79	81	83	86	88	89	89	91	91	94	94	97	100	101	104	104	104	104	104	106	107	107	109	109	111	111	113	115	117	119	122	124	126	126	127	127	127	127	129	130	133	135	135	138	141	143	144	146	146	148	148	149	149	149	152	152	153	156	158	160	163	163	164	165	165	166	167	167	169	171	174	175	175	178	178	179	181	184	187	187	189	190	190	193	196	196	198	198	199	200	203	206	209	211	211	211	211	214	217	219	220	222	225	226	226	227	227	228	230	230	230	233	233	235	237	237	237	238	238	241	242	242	245	246	246	247	250	251	253	258	260	261	261	264	265	267	269	271	274	277	277	280	281	282	283	285	286	289	291	292	293	293	295	297	297	298	300	302	304	305	307	309	309	311	314	315	315	316	316	318	321	321	322	325	327	327	327	328	328	330	333	334	336	337	339	340	340	342	346	348	348	349	349	352	353	354	355	357	359	361	363	364	365	367	371	373	375	375	375	376	377	380	383	385	387	388	389	391	393	395	395	397	398
Bases:	tcagtattcgtgtcagCATGCTGCCTCCCGTAGGAGTTTGGACCGTGTCTCAGTTCCAATGTGGGGGACCTTCCTCTCAGAACCCCTATCCATCGAAGGTTTGGTGAGCCGTTACCTCACCAACTGCCTAATGGAACGCATCCCCATCGATAACCGAAATTCTTTAATAATAACACnngttgtntcattgtactatcgttnttanagtcttnaccggnttatcccggntagtcggnaggttgtactcncgtncncccgtgcncgtcnccta
Quality Scores:	32	32	32	32	32	32	32	32	32	32	32	32	32	32	32	32	32	32	32	32	32	32	40	40	40	40	32	32	32	40	40	40	40	40	40	40	40	40	40	40	40	40	40	40	40	40	40	40	40	32	32	32	32	32	32	32	32	32	32	32	32	32	31	31	31	31	31	32	32	32	32	32	32	32	32	32	32	32	32	32	32	32	21	21	21	21	32	32	32	32	32	32	32	32	32	32	32	32	32	32	32	32	32	32	32	32	32	32	28	28	28	30	32	32	32	32	32	32	32	32	32	32	32	32	32	32	32	32	32	32	32	32	32	32	32	32	32	32	32	32	32	32	32	32	32	32	32	32	32	32	32	32	32	32	32	32	17	17	17	32	32	32	30	30	28	29	29	21	24	20	15	13	13	13	13	18	0	0	13	13	13	13	13	0	13	19	19	13	13	13	13	19	23	23	23	21	18	12	12	12	0	12	12	12	0	12	12	12	21	12	12	0	12	12	12	12	12	0	12	12	11	20	17	15	15	11	11	0	11	17	11	11	11	11	11	0	11	11	11	11	11	11	15	17	17	11	11	0	11	11	11	0	11	0	11	11	11	11	11	15	11	0	11	18	18	11	0	11	11	16	18
"""

if __name__ == "__main__":
    main()
