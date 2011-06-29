#!/usr/bin/env python
# File created on 19 Mar 2011
from __future__ import division

__author__ = "Justin Kucyznski"
__copyright__ = "Copyright 2011, The QIIME Project"
__credits__ = ["Justin Kucyznski"]
__license__ = "GPL"
__version__ = "1.3.0"
__maintainer__ = "Justin Kucyznski"
__email__ = "justinak@gmail.com"
__status__ = "Release"

from cogent.util.unit_test import TestCase, main
from cogent.parse.tree import DndParser

from qiime.parse import parse_otu_table
from qiime.util import get_qiime_scripts_dir, load_qiime_config

import tempfile
import string
import random
import os
import shutil
import subprocess
import numpy

import qiime.simsam

class SimsamTests(TestCase):

    def setUp(self):
        self.dirs_to_remove = []

    def tearDown(self):
        for d in self.dirs_to_remove:
            if os.path.exists(d):
                shutil.rmtree(d)

    def test_script(self):
        """ test the whole simsam script
        """
        qiime_config = load_qiime_config()
        tempdir = qiime_config['temp_dir'] or tempfile.gettempdir()
        maindir = os.path.join(tempdir,
         ''.join(random.choice(string.ascii_letters + string.digits) \
         for x in range(10)))


        os.makedirs(maindir)
        self.dirs_to_remove.append(maindir)
        otuf = os.path.join(maindir,'otuf')
        treef = os.path.join(maindir,'treef')

        otufh = open(otuf,'w')
        otufh.write(tutorial_otu_table)
        otufh.close()

        treefh = open(treef,'w')
        treefh.write(tutorial_tree)
        treefh.close()

        scripts_dir = get_qiime_scripts_dir()
        cmd = scripts_dir+'/simsam.py -i %s -t %s -o %s -d .003 -n 3  ' %\
         (otuf, treef, maindir+'/simsam_out.txt')
        proc = subprocess.Popen(cmd,shell=True, 
            stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        scriptout, scripterr = proc.communicate()

        if scriptout:
            raise RuntimeError('script returned stdout: ' + scriptout)
        if scripterr:
            raise RuntimeError('script returned stderr: ' + scripterr)

        num_replicates = 3 # ensure this matches cmd above

        res = open(maindir+'/simsam_out.txt','U')
        sample_ids, otu_ids, otu_mtx, lineages = parse_otu_table(open(otuf,'U'))
        res_sample_ids, res_otu_ids, res_otu_mtx, res_lineages =\
            parse_otu_table(res)

        # 3 samples per input sample
        self.assertEqual(num_replicates*len(sample_ids),len(res_sample_ids))

        # sample_ids have correct naming and order
        for i in range(len(sample_ids)):
            for j in range(num_replicates):
                exp = sample_ids[i] +'.'+str(j)
                self.assertEqual(res_sample_ids[i*num_replicates+j],exp)

        # same total sequences in each replicate sample
        num_output_samples = len(sample_ids) * num_replicates
        for i in range(num_replicates):
            self.assertEqual(res_otu_mtx[:,
             range(i,num_output_samples,num_replicates)].sum(0),
             otu_mtx.sum(0))

        # would be nice to test that result otu table doesn't match input,
        # but not sure how probable that is, and don't want stochastic failures

    def test_script_nochange(self):
        """ simsam script with 0 distance should just replicate input samples
        """
        qiime_config = load_qiime_config()
        tempdir = qiime_config['temp_dir'] or tempfile.gettempdir()
        maindir = os.path.join(tempdir,
         ''.join(random.choice(string.ascii_letters + string.digits) \
         for x in range(10)))


        os.makedirs(maindir)
        self.dirs_to_remove.append(maindir)
        otuf = os.path.join(maindir,'otuf')
        treef = os.path.join(maindir,'treef')

        otufh = open(otuf,'w')
        otufh.write(tutorial_otu_table)
        otufh.close()

        treefh = open(treef,'w')
        treefh.write(tutorial_tree)
        treefh.close()

        scripts_dir = get_qiime_scripts_dir()
        cmd = scripts_dir+'/simsam.py -i %s -t %s -o %s -d 0 -n 3  ' %\
         (otuf, treef, maindir+'/simsam_out.txt')
        proc = subprocess.Popen(cmd,shell=True, 
            stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        scriptout, scripterr = proc.communicate()

        if scriptout:
            raise RuntimeError('script returned stdout: ' + scriptout)
        if scripterr:
            raise RuntimeError('script returned stderr: ' + scripterr)

        num_replicates = 3 # ensure this matches cmd above
        res = open(maindir+'/simsam_out.txt','U')
        sample_ids, otu_ids, otu_mtx, lineages = parse_otu_table(open(otuf,'U'))
        res_sample_ids, res_otu_ids, res_otu_mtx, res_lineages =\
            parse_otu_table(res)

        # 3 samples per input sample
        self.assertEqual(num_replicates*len(sample_ids),len(res_sample_ids))

        # sample_ids have correct naming and order
        for i in range(len(sample_ids)):
            for j in range(num_replicates):
                exp = sample_ids[i] +'.'+str(j)
                self.assertEqual(res_sample_ids[i*num_replicates+j],exp)

        # same otu ids
        self.assertEqual(res_otu_ids, otu_ids)

        # same otu table, just replicated thrice
        # note this requires the same sorting of otus, input is correct sorting
        num_output_samples = len(sample_ids) * num_replicates
        for i in range(num_replicates):
            self.assertEqual(res_otu_mtx[:,
             range(i,num_output_samples,num_replicates)],
             otu_mtx)


    def test_sim_otu_table(self):
        """ simulated otu table should be right order, number of seqs

        tree looks like:
                  /-A
                 |
        ---------|--B
                 |
                 |          /-C
                  \--------|
                            \-D
        """
        sample_ids = ['samB','samA']
        otu_ids = ['C','A']
        otu_mtx = numpy.array([ [3,9],
                                [5,0],
                                ])
        otu_metadata = ['otu_C is cool','']
        tree = DndParser("(A:0.1,B:0.2,(C:0.3,D:0.4):0.5);")
        num_replicates = 3
        dissimilarity = 0.15
        res_sam_names, res_otus, res_otu_mtx, res_otu_metadata = \
         qiime.simsam.sim_otu_table(sample_ids, otu_ids, otu_mtx, otu_metadata,
         tree, num_replicates, dissimilarity)
        
        # dissim is too small to change otu C, it should always be there
        # with at least original # seqs, maybe more
        c_index = res_otus.index('C')
        c_row = res_otu_mtx[c_index]
        self.assertEqual(res_otu_metadata[c_index],'otu_C is cool')
        self.assertGreaterThan(c_row,[2,2,2,8,8,8])

        # order of samples should remain the same as input,
        # and eash replicate sample
        # should have same number or sequences as input
        for i in range(len(sample_ids)):
            for j in range(num_replicates):
                self.assertEqual(otu_mtx[:,i].sum(),
                 res_otu_mtx[:,num_replicates*i+j].sum())

    def test_get_new_otu_id_small(self):
        """ small dissim should return old tip id"""
        tree = DndParser("(A:0.1,B:0.2,(C:0.3,D:0.4):0.5);")
        res = qiime.simsam.get_new_otu_id(old_otu_id='A', tree=tree, dissim=.05)
        self.assertEqual(res,'A')

    def test_get_new_otu_id_large(self):
        """  w/ large dissim, should at least sometimes return other tip"""
        tree = DndParser("(A:0.1,B:0.2,(C:0.3,D:0.4):0.5);")
        results = []
        for i in range(1000):
            results.append(qiime.simsam.get_new_otu_id(old_otu_id='D', 
            tree=tree, dissim=.6))
        self.assertContains(results,'C')
        self.assertContains(results,'D')
        self.assertNotContains(results,'A')
        self.assertNotContains(results,'B')

    def test_combine_sample_dicts(self):
        """ combining sample dicts should give correct otu table and sorting
        """
        d1 = {'otu2':0,'otu1':3}
        d2 = {'otu4':5}
        d3 = {}
        res_otu_mtx, res_otu_ids = qiime.simsam.combine_sample_dicts([d1,d2,d3])
        exp_otu_ids = ['otu1','otu2','otu4']
        exp_otu_mtx = numpy.array([ [3,0,0],
                                    [0,0,0],
                                    [0,5,0],
                                    ])
        self.assertEqual(res_otu_ids, exp_otu_ids)
        self.assertEqual(res_otu_mtx, exp_otu_mtx)


tutorial_tree = """(((((381:0.0213,(214:0.03728,253:0.00015)0.945:0.03224)0.763:0.00483,((269:0.02693,(231:0.00509,(105:0.01425,141:0.02641)0.846:0.01405)0.428:0.00519)0.622:0.00014,404:0.00524)0.795:0.00514)0.773:0.00508,(131:0.00518,(33:0.01631,((284:0.00828,(176:0.03098,388:0.01236)0.901:0.02175)0.885:0.01273,52:0.01046)0.743:0.00498)0.924:0.01603)0.779:0.00511)0.772:0.00014,153:0.00507)0.753:0.00602,(223:0.03237,(172:0.01733,81:0.03834)0.224:0.00414)0.845:0.01076,(136:0.00627,((((265:0.01557,200:0.00517)0.674:0.00014,((204:0.00015,(339:0.01613,(322:0.01633,268:0.01643)0.569:0.0107)0.885:0.00016)0.840:0.00527,((((((((280:0.02348,(395:0.00015,(48:0.03014,((30:0.02665,316:0.01921)0.813:0.01152,215:0.0242)0.850:0.01191)0.320:0.00016)0.912:0.02431)0.694:0.01482,(115:0.01526,364:0.08211)0.879:0.03637)0.677:0.03567,((((((162:0.06933,59:0.02113)0.991:0.08563,(308:0.02061,43:0.03488)0.894:0.04949)0.911:0.05006,(((344:0.00015,(146:0.00015,377:0.01634)0.924:0.0108)0.918:0.01069,((201:0.011,240:0.04792)1.000:0.00015,(61:0.00015,96:0.00523)0.781:0.00514)0.828:0.01056)0.809:0.00016,196:0.04505)0.213:0.00014)0.650:0.00529,(((161:0.01191,(390:0.04307,37:0.03893)0.933:0.03396)0.814:0.01401,68:0.04946)0.953:0.03303,((341:0.01127,393:0.02765)0.941:0.02238,(82:0.01112,(350:0.01141,(156:0.01636,356:0.00015)0.863:0.02214)0.946:0.02475)0.748:0.00565)0.761:0.00968)0.748:0.00836)0.927:0.0224,271:0.05902)0.753:0.00511,(((((((217:0.03796,379:0.00016)0.973:0.05805,(299:0.08963,(382:0.06426,((317:0.00016,((205:0.00532,264:0.03867)0.939:0.01605,(194:0.03374,(32:0.01052,(348:0.02212,157:0.02743)1.000:0.00014)0.868:0.02793)0.745:0.00531)0.336:0.01061)0.789:0.00604,334:0.02104)0.598:0.01527)0.687:0.00354)0.836:0.01564)0.811:0.01617,(292:0.06237,84:0.02159)0.934:0.04776)0.864:0.02103,301:0.06716)0.698:0.0046,272:0.00539)0.809:0.0115,88:0.05965)0.860:0.01208,(276:0.01065,279:0.03443)0.891:0.01124)0.090:0.00014)0.924:0.03938)0.953:0.05227,281:0.02828)0.691:0.00622,25:0.01213)0.727:0.00397,((261:0.01613,((147:0.01555,20:0.00016)0.967:0.02125,(107:0.01089,349:0.03426)0.757:0.00478)0.750:0.00518)0.799:0.0052,(259:0.01616,63:0.01053)0.764:0.00523)0.792:0.00511)1.000:0.00016,(72:0.05949,(1:0.01425,67:0.0377)0.751:0.00762)0.867:0.01609)0.807:0.00507,((49:0.01645,116:0.01633)0.736:0.00514,(398:0.00515,(((180:0.04458,99:0.0328)0.913:0.02521,(((410:0.05589,(((150:0.04425,(170:0.03163,((250:0.00693,331:0.00435)1.000:0.10845,357:0.01319)0.850:0.0225)0.879:0.02887)0.749:0.00795,(((((23:0.00919,248:0.08024)0.405:0.03691,(358:0.05635,369:0.07223)0.978:0.09469)0.888:0.05975,(234:0.07249,8:0.00016)0.712:0.01829)0.976:0.07916,(((275:0.094,(((114:0.0269,302:0.02202)0.985:0.06964,(213:0.06889,42:0.03436)0.415:0.01928)0.795:0.02064,((110:0.05188,342:0.01457)0.967:0.08524,((123:0.02756,343:0.0481)0.800:0.01738,((298:0.03283,(124:0.02507,6:0.03351)0.781:0.01076)0.939:0.03194,309:0.04124)0.820:0.01321)0.985:0.0961)0.928:0.06559)0.902:0.03886)0.684:0.03217,373:0.06838)0.909:0.03592,((290:0.02673,380:0.00015)1.000:0.16099,(((90:0.09952,192:0.10171)0.679:0.01316,(326:0.03972,45:0.09053)0.965:0.05309)0.115:0.00014,(375:0.00015,(221:0.00071,278:0.05255)1.000:0.08313)1.000:0.10921)0.623:0.0222)0.892:0.03509)0.465:0.00015)0.980:0.05443,(((306:0.08813,385:0.14214)0.269:0.00862,((256:0.01776,(273:0.07543,69:0.01333)0.591:0.02343)0.883:0.02549,((132:0.02365,219:0.01597)0.897:0.02388,(100:0.01243,50:0.0237)0.226:0.01766)0.961:0.04348)0.848:0.01577)0.998:0.08323,(241:0.23207,(130:0.24778,(53:0.12887,(129:0.07692,318:0.01288)0.900:0.04845)0.817:0.02143)0.888:0.05464)0.657:0.01537)0.822:0.01876)0.828:0.01549)0.773:0.01019,((98:0.12681,((148:0.0294,391:0.00571)0.989:0.07803,(389:0.10107,(252:0.00014,362:0.01104)0.964:0.06682)0.834:0.03217)0.762:0.0152)0.524:0.0181,(0:0.0483,(135:0.01151,(300:0.0175,(274:0.04561,((((166:0.02935,355:0.00015)0.833:0.00565,41:0.00014)0.807:0.00586,(226:0.01038,92:0.0044)0.792:0.00425)0.961:0.03236,((360:0.01752,7:0.0182)0.748:0.00495,(368:0.02316,288:0.01783)0.759:0.00622)0.707:0.00573)0.841:0.00015)0.949:0.02275)0.745:0.00559)0.855:0.02344)0.876:0.03532)0.885:0.02567)0.752:0.00645)0.782:0.00969,(((((((178:0.01576,(230:0.02704,64:0.02146)0.869:0.0108)0.809:0.01014,((122:0.00448,354:0.0225)0.855:0.01127,(333:0.01086,406:0.01648)0.748:0.00433)0.789:0.00624)0.171:0.00516,((416:0.04298,(400:0.01045,74:0.01051)0.923:0.00014)0.862:0.02166,(307:0.04097,(260:0.03574,335:0.0434)0.747:0.00875)0.916:0.02837)0.843:0.00987)0.804:0.00016,((237:0.09447,((370:0.01631,(319:0.04803,(60:0.01986,405:0.01742)0.560:0.01574)0.898:0.01971)0.918:0.01584,(384:0.02116,(245:0.01047,(177:0.0051,(183:0.03226,413:0.00014)0.826:0.00518)0.777:0.00501)0.923:0.0158)0.622:0.00016)0.685:0.00099)0.224:0.02406,((22:0.03142,5:0.06696)0.870:0.03448,47:0.0347)0.763:0.01052)0.847:0.01209)0.743:0.00534,((((62:0.00137,(121:0.00016,78:0.04376)1.000:0.10609)0.942:0.0378,(311:0.05626,407:0.06902)0.944:0.04614)0.703:0.00608,(((188:0.01993,202:0.02611)0.914:0.02118,(328:0.0273,337:0.00015)0.815:0.01019)0.852:0.01169,(330:0.03441,((386:0.13035,(392:0.00544,(321:0.02191,4:0.01061)0.763:0.0052)0.932:0.00014)0.671:0.01096,145:0.01556)0.829:0.01073)0.735:0.00529)0.840:0.01052)0.849:0.01531,(262:0.0683,((310:0.05551,((83:0.01296,(127:0.01909,212:0.01393)0.090:0.00499)0.876:0.01352,(104:0.00014,171:0.01061)0.895:0.01717)0.877:0.02683)0.940:0.03929,(119:0.0152,179:0.00197)0.889:0.02843)0.066:0.01551)0.839:0.01374)0.820:0.01069)0.869:0.01061,(((293:0.01741,168:0.04514)0.046:0.01491,345:0.03334)0.248:0.01629,(31:0.04727,97:0.04999)0.915:0.03556)0.811:0.01631)0.010:0.00016,(((94:0.0671,(108:0.00014,229:0.06991)0.630:0.01827)0.982:0.06031,(143:0.02201,((((((198:0.02745,(140:0.14724,75:0.02831)0.817:0.0209)0.851:0.01902,(((282:0.06783,54:0.00015)0.952:0.03641,((313:0.03746,80:0.00524)0.872:0.0215,2:0.07468)0.862:0.02589)0.916:0.03914,((367:0.0099,(((128:0.0425,((111:0.06727,11:0.00495)0.974:0.02953,283:0.02606)0.504:0.00357)0.862:0.02044,(289:0.04546,(399:0.00319,((((152:0.00014,19:0.06307)0.992:0.03752,154:0.00016)0.786:0.00014,134:0.06945)0.997:0.06109,51:0.00014)0.994:0.04556)0.353:0.00583)0.482:0.00828)0.933:0.03536,112:0.07957)0.734:0.00733)0.962:0.08492,403:0.10375)0.869:0.0525)0.894:0.03949)0.645:0.00925,((((287:0.00534,15:0.05518)0.920:0.03189,(((304:0.00508,409:0.00015)0.991:0.00014,(120:0.00015,(57:0.04309,56:0.0156)0.759:0.00015)0.902:0.01019)0.339:0.01644,173:0.094)0.787:0.01131)1.000:0.07731,(236:0.00625,((26:0.04569,(((351:0.005,(27:0.03624,(137:0.01569,(314:0.00015,408:0.03277)0.991:0.03257)0.806:0.00498)0.851:0.00588)0.928:0.01791,((133:0.04374,(227:0.00527,(412:0.00014,(175:0.00507,((95:0.01566,210:0.00014)0.438:0.01045,191:0.00016)0.781:0.00518)0.815:0.00508)0.859:0.01021)0.745:0.00667)0.735:0.01956,((((12:0.01588,415:0.01701)0.121:0.03139,(73:0.04886,(17:0.00016,(46:0.02083,378:0.01021)0.886:0.01027)0.785:0.019)0.719:0.02118)0.774:0.01959,329:0.01522)0.777:0.01121,(((286:0.00722,(394:0.01596,(372:0.00015,225:0.0446)0.884:0.0109)0.929:0.02558)0.584:0.00985,218:0.02283)0.888:0.01478,159:0.02121)0.739:0.00866)0.851:0.01129)0.728:0.00602)0.866:0.01998,93:0.04869)0.604:0.00297)0.648:0.01633,199:0.06704)0.788:0.01956)0.371:0.01052)0.827:0.01491,((244:0.0262,(126:0.00015,163:0.03192)0.984:0.04377)0.817:0.01306,((216:0.00014,(86:0.02257,(21:0.01127,34:0.01066)0.859:0.01088)0.622:0.017)0.998:0.19434,(233:0.00244,(182:0.01898,(239:0.02877,267:0.00015)0.839:0.01438)0.999:0.09419)0.975:0.15234)0.877:0.07457)0.893:0.0244)0.821:0.02013)0.998:0.10422,(195:0.10508,((249:0.0368,(336:0.04596,((263:0.02407,(277:0.01295,190:0.03788)0.823:0.01671)0.698:0.0068,197:0.01756)0.309:0.01631)0.860:0.01866)0.926:0.02656,(303:0.04293,(113:0.04423,347:0.04295)0.930:0.03972)0.885:0.02484)0.701:0.00015)0.902:0.03629)0.841:0.02905,(246:0.00014,(125:0.03009,184:0.0229)0.998:0.07478)0.999:0.10301)0.936:0.04978,((247:0.04204,((((((238:0.01393,(109:0.01081,39:0.02762)0.769:0.00519)0.758:0.00702,(257:0.01539,85:0.07408)0.746:0.00558)0.755:0.01039,(363:0.04294,155:0.00015)0.943:0.02426)0.894:0.01745,266:0.00586)0.948:0.03346,55:0.02705)0.739:0.00453,203:0.00015)0.855:0.02077)0.995:0.07638,327:0.00745)0.921:0.03692)0.553:0.01549)0.970:0.05544)0.858:0.02855,338:0.08163)0.892:0.03304)0.759:0.00673)0.945:0.02495,((((((((102:0.04317,36:0.02415)0.964:0.03758,65:0.00505)0.822:0.01078,366:0.00016)0.811:0.01537,(315:0.01071,((151:0.0,160:0.0):0.00016,340:0.00014)0.842:0.01037)0.951:0.02712)0.724:0.00057,(185:0.04527,(207:0.01304,76:0.00341)0.949:0.03474)0.845:0.0196)0.871:0.0106,(371:0.02805,(164:0.0104,242:0.02179)0.758:0.0052)0.771:0.00538)0.841:0.01097,174:0.13953)0.831:0.01033,(144:0.01866,(3:0.01578,312:0.00015)0.785:0.00532)0.780:0.00615)0.752:0.00572)0.667:0.00244)0.268:0.00339,((101:0.04199,77:0.00334)0.965:0.0345,((14:0.01106,294:0.00502)0.891:0.01811,(285:0.01062,397:0.01076)0.758:0.00896)0.163:0.01034)0.850:0.01331)0.563:0.00537)0.800:0.00519)0.930:0.00016)0.759:0.01023)1.000:0.00014)0.850:0.00015,(243:0.03373,220:0.01032)0.888:0.011)0.540:0.00014,(189:0.02629,(((139:0.0155,186:0.01757)0.420:0.01444,(((((((165:0.0059,58:0.03297)0.779:0.02132,((222:0.01678,(323:0.02243,44:0.04081)0.819:0.01102)0.063:0.00015,(106:0.03989,149:0.02047)0.775:0.01298)0.706:0.0074)0.957:0.03281,((((258:0.04247,87:0.0123)0.500:0.01067,235:0.00735)0.645:0.00296,208:0.00505)1.000:0.00015,((18:0.00454,(((10:0.04233,(414:0.00016,(142:0.01127,66:0.03479)0.756:0.00498)0.726:0.00685)0.486:0.01639,181:0.00014)0.784:0.00501,(167:0.01463,(320:0.00885,402:0.00881)0.791:0.00014)0.839:0.01499)0.773:0.00524)0.893:0.01079,(169:0.00517,(295:0.01586,297:0.03792)0.262:0.00016)0.778:0.00521)0.818:0.00528)0.764:0.01062)0.767:0.00486,70:0.00512)0.766:0.00495,(((332:0.00016,((325:0.01591,(383:0.00014,(361:0.01642,(138:0.04133,(158:0.0036,224:0.00657)0.840:0.01972)0.769:0.00881)0.777:0.00496)0.882:0.01036)0.752:0.00492,(24:0.03974,((((254:0.00541,(251:0.00015,(324:0.02187,((117:0.0052,(374:0.03165,270:0.02362)0.731:0.00708)0.791:0.00525,13:0.01621)0.757:0.00511)0.607:0.01283)0.889:0.0192)0.852:0.01583,305:0.01647)0.948:0.00015,211:0.00015)0.419:0.00016,(103:0.01686,209:0.05269)0.861:0.01595)0.937:0.01635)0.756:0.00523)0.878:0.01048)0.776:0.00238,(365:0.03251,((38:0.04434,79:0.00014)0.758:0.00016,(296:0.043,9:0.00518)0.693:0.0162)0.508:0.00805)0.766:0.00767)0.764:0.00313,(((359:0.02181,(16:0.04469,(232:0.01621,(118:0.03421,(29:0.01612,353:0.01494)0.293:0.01034)0.864:0.01326)0.747:0.01394)0.724:0.0072)0.911:0.01681,387:0.02755)0.761:0.00523,(346:0.01957,(376:0.04072,71:0.0547)0.829:0.0181)0.750:0.00673)0.823:0.01037)0.774:0.0054)0.789:0.005,(((228:0.00529,((401:0.02214,((187:0.00532,411:0.00526)0.801:0.00583,((89:0.027,193:0.00014)0.787:0.00524,91:0.01618)0.743:0.0045)0.548:0.00548)0.825:0.016,40:0.02807)0.778:0.00992)0.824:0.01011,255:0.05012)0.966:0.00014,(352:0.01585,396:0.00014)0.784:0.02134)0.880:0.0107)0.901:0.0194,(35:0.0209,(206:0.00836,291:0.06414)0.439:0.00793)0.753:0.00846)0.763:0.00968)0.942:0.02851,28:0.0208)0.742:0.01057)0.781:0.00811)0.802:0.02029)0.750:0.01578);"""

tutorial_otu_table = """#Full OTU Counts
#OTU ID	PC.354	PC.355	PC.356	PC.481	PC.593	PC.607	PC.634	PC.635	PC.636	Consensus Lineage
0	0	0	0	0	0	0	0	1	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
1	0	0	0	0	0	1	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
2	0	0	0	0	0	0	0	0	1	Root;Bacteria;Bacteroidetes;Bacteroidetes;Bacteroidales;Porphyromonadaceae;Parabacteroides
3	2	1	0	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae";"Lachnospiraceae Incertae Sedis"
4	1	0	0	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
5	0	0	0	0	0	0	0	0	1	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
6	0	0	0	0	0	0	0	1	0	Root;Bacteria;Actinobacteria;Actinobacteria
7	0	0	2	0	0	0	0	0	2	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales
8	1	1	0	2	4	0	0	0	0	Root;Bacteria;Firmicutes;"Bacilli";"Lactobacillales";Lactobacillaceae;Lactobacillus
9	0	0	2	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
10	0	1	0	0	0	0	0	0	0	Root;Bacteria
11	0	0	0	0	0	0	1	0	0	Root;Bacteria;Bacteroidetes;Bacteroidetes;Bacteroidales;Bacteroidaceae;Bacteroides
12	0	0	0	0	0	0	1	0	0	Root;Bacteria;Bacteroidetes
13	1	0	0	1	0	1	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales
14	0	0	1	1	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
15	0	0	0	0	1	0	0	0	0	Root;Bacteria
16	1	0	2	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
17	0	0	0	1	0	0	4	10	37	Root;Bacteria;Bacteroidetes
18	0	1	0	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
19	0	0	0	0	0	0	0	0	1	Root;Bacteria;Bacteroidetes
20	0	0	0	0	1	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
21	0	0	0	0	0	0	2	3	2	Root;Bacteria;Bacteroidetes
22	0	0	0	0	2	0	1	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
23	14	1	14	1	0	0	0	0	0	Root;Bacteria;Firmicutes;"Bacilli";"Lactobacillales";Lactobacillaceae;Lactobacillus
24	1	0	0	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
25	0	0	0	1	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae";"Lachnospiraceae Incertae Sedis"
26	0	0	0	0	0	0	0	1	1	Root;Bacteria;Bacteroidetes
27	0	0	0	0	0	0	0	0	1	Root;Bacteria;Bacteroidetes
28	0	1	0	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
29	6	0	4	0	2	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
30	0	0	0	0	0	1	0	0	0	Root;Bacteria;Firmicutes
31	1	0	0	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales
32	0	0	0	0	1	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Ruminococcaceae"
33	0	0	0	1	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
34	0	0	0	0	0	0	8	10	2	Root;Bacteria
35	1	0	1	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
36	1	0	1	0	0	0	0	1	1	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
37	0	0	0	0	0	1	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales
38	0	0	1	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales
39	0	0	0	0	0	0	0	1	0	Root;Bacteria;Bacteroidetes;Bacteroidetes;Bacteroidales;Rikenellaceae;Alistipes
40	0	0	1	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
41	0	0	1	0	0	0	0	1	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales
42	0	0	0	0	0	1	0	0	0	Root;Bacteria;Firmicutes
43	0	0	0	0	0	1	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales
44	0	0	1	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales
45	1	0	0	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Erysipelotrichi";"Erysipelotrichales";Erysipelotrichaceae;Coprobacillus
46	0	0	0	0	0	0	0	0	1	Root;Bacteria;Bacteroidetes
47	0	0	0	1	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
48	0	0	0	0	1	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales
49	0	0	0	1	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
50	0	1	0	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
51	0	1	0	0	0	0	0	0	0	Root;Bacteria;Bacteroidetes;Bacteroidetes;Bacteroidales;Bacteroidaceae;Bacteroides
52	0	2	0	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
53	0	0	0	0	0	0	2	0	1	Root;Bacteria;Proteobacteria;Deltaproteobacteria
54	0	0	0	0	0	0	5	0	0	Root;Bacteria;Bacteroidetes;Bacteroidetes;Bacteroidales;Porphyromonadaceae;Parabacteroides
55	0	0	0	0	0	0	1	0	0	Root;Bacteria;Bacteroidetes;Bacteroidetes;Bacteroidales;Rikenellaceae;Alistipes
56	0	0	0	0	0	1	0	0	0	Root;Bacteria;Bacteroidetes
57	0	0	0	0	0	0	0	1	0	Root;Bacteria
58	1	0	1	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
59	0	0	0	0	0	0	0	0	1	Root;Bacteria;Deferribacteres;Deferribacteres;Deferribacterales;Deferribacteraceae;Mucispirillum
60	0	0	0	0	0	0	0	1	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
61	0	0	1	0	0	0	0	1	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Ruminococcaceae"
62	0	0	1	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales
63	1	0	1	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
64	0	0	0	0	0	0	0	0	1	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales
65	0	0	0	6	0	0	0	1	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
66	0	0	1	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
67	0	0	1	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
68	1	0	0	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales
69	0	0	1	0	0	0	0	0	0	Root;Bacteria
70	0	0	0	0	0	1	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
71	0	0	1	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
72	0	0	0	0	0	1	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
73	0	0	0	0	0	5	0	0	0	Root;Bacteria;Bacteroidetes
74	0	0	0	1	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
75	1	0	1	0	0	0	0	0	0	Root;Bacteria;Bacteroidetes
76	0	0	0	1	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
77	0	0	0	1	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
78	1	0	1	1	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales
79	2	3	8	0	1	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
80	0	0	0	0	0	0	0	0	1	Root;Bacteria;Bacteroidetes;Bacteroidetes;Bacteroidales;Porphyromonadaceae;Parabacteroides
81	1	0	0	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae";"Lachnospiraceae Incertae Sedis"
82	0	0	0	0	0	2	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales
83	0	0	0	1	0	0	0	1	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
84	1	0	0	0	0	0	0	2	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Ruminococcaceae";Ruminococcus
85	0	0	0	0	0	0	0	0	1	Root;Bacteria;Bacteroidetes;Bacteroidetes;Bacteroidales;Rikenellaceae;Alistipes
86	0	0	0	0	0	0	0	1	0	Root;Bacteria
87	0	0	1	0	0	2	0	1	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
88	0	0	0	0	0	0	0	1	0	Root;Bacteria
89	0	0	1	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
90	0	0	0	9	0	0	3	0	0	Root;Bacteria;Firmicutes;"Erysipelotrichi";"Erysipelotrichales";Erysipelotrichaceae;Turicibacter
91	0	0	0	1	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae";Butyrivibrio
92	0	0	0	0	0	0	1	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales
93	0	0	0	0	0	0	2	1	0	Root;Bacteria;Bacteroidetes
94	0	0	0	0	0	0	0	1	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales
95	0	0	0	2	0	0	0	0	0	Root;Bacteria;Bacteroidetes
96	0	0	0	1	0	1	0	1	1	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Ruminococcaceae"
97	0	0	0	0	0	1	0	0	0	Root;Bacteria
98	0	0	0	0	0	0	0	1	0	Root;Bacteria
99	0	0	0	1	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
100	0	0	0	1	0	0	0	0	0	Root;Bacteria
101	0	0	0	3	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales
102	0	1	0	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
103	0	1	0	0	0	0	1	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
104	0	0	0	0	0	1	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales
105	0	1	0	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
106	0	0	0	0	0	1	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales
107	0	0	0	0	0	1	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
108	0	0	0	0	0	0	1	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;Incertae Sedis XIII;Anaerovorax
109	0	0	0	1	0	0	1	5	2	Root;Bacteria;Bacteroidetes;Bacteroidetes;Bacteroidales;Rikenellaceae;Alistipes
110	0	0	0	0	0	2	0	0	0	Root;Bacteria;Actinobacteria;Actinobacteria;Coriobacteridae;Coriobacteriales;Coriobacterineae;Coriobacteriaceae;Olsenella
111	0	0	0	0	0	0	1	0	0	Root;Bacteria;Bacteroidetes;Bacteroidetes;Bacteroidales;Bacteroidaceae;Bacteroides
112	0	0	0	0	0	0	1	0	0	Root;Bacteria;Bacteroidetes;Bacteroidetes;Bacteroidales;Bacteroidaceae;Bacteroides
113	0	0	0	0	0	1	0	0	0	Root;Bacteria
114	0	0	0	0	0	1	0	0	0	Root;Bacteria
115	0	0	0	0	0	1	0	0	0	Root;Bacteria;Firmicutes
116	0	1	0	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae";"Lachnospiraceae Incertae Sedis"
117	1	0	2	0	0	6	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales
118	0	0	0	1	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
119	0	0	0	0	0	0	0	1	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales
120	1	3	1	2	1	9	2	4	5	Root;Bacteria;Bacteroidetes
121	0	0	0	0	0	0	0	1	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales
122	0	0	0	1	0	2	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
123	0	0	0	0	0	0	1	0	0	Root;Bacteria;Actinobacteria;Actinobacteria;Coriobacteridae;Coriobacteriales;Coriobacterineae;Coriobacteriaceae
124	0	0	0	0	0	0	1	0	0	Root;Bacteria;Actinobacteria;Actinobacteria;Coriobacteridae;Coriobacteriales;Coriobacterineae;Coriobacteriaceae
125	0	0	0	0	0	0	1	0	0	Root;Bacteria;Bacteroidetes
126	0	0	2	0	0	0	0	1	0	Root;Bacteria
127	0	0	0	0	0	1	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
128	0	0	0	0	0	0	1	0	0	Root;Bacteria;Bacteroidetes;Bacteroidetes;Bacteroidales;Bacteroidaceae;Bacteroides
129	0	0	0	1	0	0	0	0	0	Root;Bacteria
130	0	0	0	0	5	2	0	0	0	Root;Bacteria;Proteobacteria;Epsilonproteobacteria;Campylobacterales;Helicobacteraceae;Helicobacter
131	0	0	1	3	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae";"Lachnospiraceae Incertae Sedis"
132	0	0	0	0	1	0	0	0	0	Root;Bacteria
133	0	0	1	0	0	0	0	0	0	Root;Bacteria
134	0	0	0	0	0	0	0	0	1	Root;Bacteria;Bacteroidetes;Bacteroidetes;Bacteroidales;Bacteroidaceae;Bacteroides
135	0	0	1	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales
136	1	0	0	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae";"Lachnospiraceae Incertae Sedis"
137	0	0	0	0	0	0	0	1	0	Root;Bacteria
138	0	0	1	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
139	1	0	0	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
140	0	0	0	0	0	0	1	3	0	Root;Bacteria
141	0	0	0	0	1	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
142	0	0	0	0	1	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
143	0	0	1	0	0	0	0	0	0	Root;Bacteria
144	0	0	0	0	0	1	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
145	0	0	2	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales
146	1	0	0	0	2	0	2	0	3	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Ruminococcaceae"
147	0	1	0	1	1	0	0	0	3	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
148	0	0	0	0	0	1	0	0	0	Root;Bacteria;Firmicutes
149	0	0	0	0	0	0	0	1	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
150	0	0	0	0	1	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales
151	0	0	0	1	0	0	0	1	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
152	0	0	0	1	0	0	1	2	19	Root;Bacteria;Bacteroidetes;Bacteroidetes;Bacteroidales;Bacteroidaceae;Bacteroides
153	0	2	1	2	0	0	1	1	1	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae";"Lachnospiraceae Incertae Sedis"
154	2	18	0	1	0	0	21	4	4	Root;Bacteria;Bacteroidetes;Bacteroidetes;Bacteroidales;Bacteroidaceae;Bacteroides
155	0	0	0	0	0	5	9	5	3	Root;Bacteria;Bacteroidetes;Bacteroidetes;Bacteroidales;Rikenellaceae;Alistipes
156	0	0	1	0	0	0	0	1	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales
157	0	0	1	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Ruminococcaceae"
158	1	0	1	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales
159	0	0	0	0	0	0	0	1	1	Root;Bacteria;Bacteroidetes
160	0	0	0	0	0	0	1	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
161	0	0	1	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales
162	0	0	0	0	0	3	5	2	6	Root;Bacteria;Deferribacteres;Deferribacteres;Deferribacterales;Deferribacteraceae;Mucispirillum
163	0	0	0	0	0	0	0	0	1	Root;Bacteria
164	0	0	0	0	0	1	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
165	2	1	1	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
166	0	0	0	0	0	0	0	1	0	Root;Bacteria
167	1	0	0	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
168	0	0	0	1	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
169	0	2	0	7	0	0	0	2	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
170	0	0	0	1	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales
171	0	0	0	1	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales
172	1	0	0	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae";"Lachnospiraceae Incertae Sedis"
173	0	0	0	0	0	1	0	0	0	Root;Bacteria
174	1	0	0	0	10	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Peptostreptococcaceae";"Peptostreptococcaceae Incertae Sedis"
175	0	0	0	0	1	0	0	0	0	Root;Bacteria;Bacteroidetes
176	0	0	0	0	0	1	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
177	0	0	0	1	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia"
178	0	0	0	2	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
179	0	0	0	1	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales
180	0	0	0	0	1	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
181	1	4	2	6	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
182	0	0	0	0	0	1	0	0	0	Root;Bacteria
183	0	0	0	0	0	0	1	0	0	Root;Bacteria;Firmicutes;"Clostridia"
184	0	0	0	1	0	0	3	1	0	Root;Bacteria;Bacteroidetes
185	0	0	0	0	0	0	0	0	1	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
186	0	0	1	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
187	0	1	0	0	0	0	0	0	1	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
188	0	0	0	0	0	0	0	1	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
189	0	0	0	1	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
190	0	0	0	0	0	0	0	1	0	Root;Bacteria
191	2	1	10	2	24	0	0	1	1	Root;Bacteria
192	0	0	0	0	0	1	0	0	0	Root;Bacteria;Firmicutes;"Bacilli";"Lactobacillales";Streptococcaceae;Streptococcus
193	0	0	0	0	0	1	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae";Butyrivibrio
194	0	0	2	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Ruminococcaceae";Acetanaerobacterium
195	0	0	0	0	0	1	0	0	0	Root;Bacteria
196	0	0	0	0	0	1	0	1	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales
197	0	1	0	0	0	0	0	0	0	Root;Bacteria
198	0	2	0	0	0	1	0	0	0	Root;Bacteria;Bacteroidetes;Bacteroidetes;Bacteroidales
199	0	0	0	0	0	1	1	0	0	Root;Bacteria
200	0	0	0	2	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
201	0	0	0	1	0	1	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales
202	0	0	0	0	0	0	1	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
203	0	2	2	4	0	5	1	5	0	Root;Bacteria;Bacteroidetes;Bacteroidetes;Bacteroidales;Rikenellaceae;Alistipes
204	1	4	0	1	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
205	0	0	0	0	0	0	0	1	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Ruminococcaceae"
206	0	1	0	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
207	0	0	0	0	0	0	0	1	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
208	0	2	0	2	0	0	0	1	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales
209	0	0	1	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales
210	0	0	0	0	0	0	0	0	1	Root;Bacteria
211	1	0	0	1	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
212	0	0	0	0	0	0	0	0	1	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
213	0	0	0	0	0	0	0	2	0	Root;Bacteria;Firmicutes
214	0	0	0	0	0	0	0	1	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
215	0	0	0	0	0	0	0	1	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales
216	0	0	0	0	0	0	0	1	0	Root;Bacteria;Bacteroidetes
217	0	0	0	0	0	2	0	1	0	Root;Bacteria
218	0	0	0	0	9	1	0	0	0	Root;Bacteria;Bacteroidetes
219	0	0	0	0	1	0	0	0	0	Root;Bacteria
220	1	0	0	0	1	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
221	0	0	0	0	0	0	0	1	0	Root;Bacteria;Firmicutes
222	0	1	0	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
223	0	0	0	0	0	0	0	2	2	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
224	0	0	0	1	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales
225	0	2	1	0	0	0	0	0	0	Root;Bacteria;Bacteroidetes
226	0	0	0	0	0	1	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales
227	0	1	2	0	9	1	1	1	3	Root;Bacteria;Bacteroidetes
228	16	0	0	0	12	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
229	0	0	0	0	0	1	1	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;Incertae Sedis XIII
230	0	0	0	1	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
231	0	19	2	0	2	0	3	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
232	0	0	0	0	0	0	1	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales
233	0	0	0	0	1	0	0	0	0	Root;Bacteria;Bacteroidetes
234	0	0	0	0	1	0	0	0	0	Root;Bacteria;Firmicutes;"Bacilli";"Lactobacillales";Lactobacillaceae;Lactobacillus
235	0	1	1	0	1	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
236	0	0	0	0	0	2	0	0	0	Root;Bacteria;Bacteroidetes;Bacteroidetes;Bacteroidales
237	0	0	0	0	1	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
238	0	0	0	0	0	0	0	1	0	Root;Bacteria;Bacteroidetes;Bacteroidetes;Bacteroidales;Rikenellaceae;Alistipes
239	0	0	0	0	0	1	0	0	0	Root;Bacteria
240	0	0	0	0	0	1	0	0	0	Root;Bacteria
241	0	0	0	0	0	0	2	0	0	Root;Bacteria;TM7;TM7_genera_incertae_sedis
242	0	0	0	0	0	0	1	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
243	0	0	0	0	0	0	1	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
244	0	0	0	0	0	0	0	0	1	Root;Bacteria;Bacteroidetes
245	0	0	0	1	0	0	0	1	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales
246	0	0	0	0	0	0	0	1	0	Root;Bacteria
247	0	0	1	0	0	0	0	0	0	Root;Bacteria;Bacteroidetes
248	1	0	0	1	0	0	0	0	0	Root;Bacteria;Firmicutes;"Bacilli";"Lactobacillales";Lactobacillaceae;Lactobacillus
249	1	0	0	0	0	0	0	0	0	Root;Bacteria
250	1	0	0	0	0	0	0	1	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales
251	0	0	0	1	4	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales
252	0	0	0	1	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Ruminococcaceae"
253	0	0	0	0	2	0	0	5	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
254	11	13	6	13	2	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
255	0	0	0	0	0	1	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
256	0	0	0	0	0	0	1	0	0	Root;Bacteria
257	0	0	0	0	0	0	5	0	0	Root;Bacteria;Bacteroidetes;Bacteroidetes;Bacteroidales;Rikenellaceae;Alistipes
258	0	0	1	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales
259	0	0	0	0	0	0	0	1	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
260	0	0	0	0	0	0	0	1	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales
261	0	0	0	0	0	0	0	1	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
262	0	1	0	0	0	0	0	0	1	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae";Bryantella
263	0	0	0	0	1	0	0	0	0	Root;Bacteria
264	0	0	0	0	0	1	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Ruminococcaceae"
265	0	0	0	0	0	2	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
266	0	0	0	2	0	0	0	0	0	Root;Bacteria;Bacteroidetes;Bacteroidetes;Bacteroidales;Rikenellaceae;Alistipes
267	1	0	0	5	17	20	0	0	0	Root;Bacteria
268	0	0	0	0	0	0	1	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
269	0	0	0	1	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae";"Lachnospiraceae Incertae Sedis"
270	0	0	1	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales
271	0	0	0	0	0	0	0	0	1	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales
272	0	0	0	1	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Ruminococcaceae"
273	0	0	0	0	0	0	1	0	0	Root;Bacteria
274	0	0	0	0	0	0	1	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales
275	0	0	0	0	0	0	1	0	0	Root;Bacteria;Verrucomicrobia;Verrucomicrobiae;Verrucomicrobiales;Verrucomicrobiaceae;Akkermansia
276	0	0	0	0	0	0	0	1	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Ruminococcaceae"
277	1	0	0	0	0	0	0	0	0	Root;Bacteria
278	0	0	0	0	0	1	0	0	0	Root;Bacteria
279	0	0	0	0	0	1	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales
280	0	1	0	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales
281	1	0	0	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae";"Lachnospiraceae Incertae Sedis"
282	0	0	0	0	0	0	2	0	0	Root;Bacteria;Bacteroidetes;Bacteroidetes;Bacteroidales;Porphyromonadaceae;Parabacteroides
283	0	0	0	0	0	0	2	1	0	Root;Bacteria;Bacteroidetes;Bacteroidetes;Bacteroidales;Bacteroidaceae;Bacteroides
284	0	0	0	1	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
285	0	0	0	0	0	0	1	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
286	0	2	3	1	4	0	5	0	4	Root;Bacteria;Bacteroidetes
287	0	0	0	0	0	0	1	1	1	Root;Bacteria;Bacteroidetes
288	0	0	0	0	0	1	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales
289	0	0	0	0	3	0	0	0	0	Root;Bacteria;Bacteroidetes;Bacteroidetes;Bacteroidales;Bacteroidaceae;Bacteroides
290	0	0	0	0	0	0	0	0	2	Root;Bacteria;Firmicutes;"Bacilli";Bacillales;"Staphylococcaceae";Staphylococcus
291	0	0	0	0	1	0	0	0	0	Root;Bacteria
292	0	0	0	0	1	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Ruminococcaceae"
293	0	0	0	0	0	1	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
294	0	1	0	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
295	29	1	10	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
296	0	0	0	0	1	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales
297	0	0	0	1	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales
298	0	0	0	0	0	0	1	0	0	Root;Bacteria;Actinobacteria;Actinobacteria
299	0	0	0	0	0	0	1	0	1	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales
300	0	0	0	0	0	1	0	0	0	Root;Bacteria;Firmicutes;"Clostridia"
301	0	0	0	0	0	0	2	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Ruminococcaceae"
302	0	0	0	0	0	1	0	0	0	Root;Bacteria
303	0	0	0	0	0	0	0	0	1	Root;Bacteria
304	0	0	0	0	0	0	0	1	0	Root;Bacteria;Bacteroidetes
305	1	0	0	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
306	0	0	0	0	0	0	0	0	1	Root;Bacteria
307	0	0	1	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
308	0	1	0	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Ruminococcaceae";"Ruminococcaceae Incertae Sedis"
309	0	0	0	1	0	0	0	0	0	Root;Bacteria;Actinobacteria;Actinobacteria;Coriobacteridae;Coriobacteriales;Coriobacterineae;Coriobacteriaceae;Denitrobacterium
310	0	0	1	0	0	0	0	0	0	Root;Bacteria
311	0	0	0	0	0	1	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales
312	0	0	1	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
313	0	1	0	0	0	0	0	0	1	Root;Bacteria;Bacteroidetes;Bacteroidetes;Bacteroidales;Porphyromonadaceae;Parabacteroides
314	0	0	1	0	0	0	0	0	0	Root;Bacteria;Bacteroidetes
315	1	3	1	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
316	0	1	0	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales
317	0	0	0	0	0	0	1	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Ruminococcaceae"
318	0	0	0	0	0	1	0	0	0	Root;Bacteria;Proteobacteria
319	0	2	1	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
320	0	0	0	1	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
321	0	0	0	0	0	0	0	0	1	Root;Bacteria
322	0	0	0	1	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
323	0	0	1	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
324	0	0	1	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
325	0	1	0	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
326	0	0	0	0	4	0	0	0	2	Root;Bacteria;Firmicutes;"Erysipelotrichi";"Erysipelotrichales";Erysipelotrichaceae;Erysipelotrichaceae Incertae Sedis
327	0	0	0	0	0	0	0	1	0	Root;Bacteria;Bacteroidetes
328	0	0	0	1	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales
329	2	2	0	1	0	0	0	0	0	Root;Bacteria;Bacteroidetes
330	0	0	1	0	0	0	0	0	0	Root;Bacteria;Firmicutes
331	0	0	0	0	1	0	0	0	0	Root;Bacteria;Firmicutes
332	0	1	0	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
333	0	0	0	0	0	6	0	3	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales
334	1	0	0	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Ruminococcaceae"
335	0	0	0	0	0	0	0	1	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
336	0	0	1	0	0	0	0	0	0	Root;Bacteria
337	0	0	0	1	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
338	0	0	0	0	0	0	0	1	0	Root;Bacteria
339	0	0	1	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
340	0	0	2	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
341	0	0	1	0	0	0	0	1	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales
342	0	0	0	0	0	1	0	0	0	Root;Bacteria
343	0	0	0	0	0	0	0	0	1	Root;Bacteria;Actinobacteria;Actinobacteria;Coriobacteridae;Coriobacteriales;Coriobacterineae;Coriobacteriaceae
344	0	0	1	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Ruminococcaceae"
345	1	0	0	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
346	0	1	0	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
347	0	0	0	1	0	0	0	0	0	Root;Bacteria
348	0	0	0	1	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Ruminococcaceae"
349	0	0	0	0	0	0	1	0	1	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
350	1	0	0	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Ruminococcaceae"
351	0	0	0	0	2	2	1	4	1	Root;Bacteria;Bacteroidetes
352	3	0	0	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
353	0	4	4	0	1	2	0	2	1	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
354	0	0	0	0	0	1	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
355	0	0	0	0	0	0	0	1	0	Root;Bacteria
356	0	0	0	0	0	1	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Ruminococcaceae"
357	0	0	0	4	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales
358	0	0	1	0	0	0	0	0	0	Root;Bacteria
359	0	0	1	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
360	0	0	1	0	0	0	0	1	1	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales
361	2	0	2	1	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
362	1	0	0	1	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Ruminococcaceae"
363	0	0	0	0	0	1	0	1	0	Root;Bacteria;Bacteroidetes;Bacteroidetes;Bacteroidales;Rikenellaceae
364	1	0	0	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales
365	0	0	0	0	0	2	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
366	0	0	0	1	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae";Roseburia
367	0	0	0	0	1	0	0	0	0	Root;Bacteria;Bacteroidetes;Bacteroidetes;Bacteroidales;Bacteroidaceae;Bacteroides
368	0	0	0	0	0	1	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales
369	0	0	0	0	0	1	0	0	0	Root;Bacteria
370	2	1	0	5	0	1	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
371	1	1	0	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
372	0	1	0	0	0	0	0	0	0	Root;Bacteria
373	0	1	0	0	0	0	3	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;Clostridiaceae;"Clostridiaceae 1";Clostridium
374	0	0	0	0	0	0	1	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
375	0	0	0	0	0	0	4	0	0	Root;Bacteria;Firmicutes;"Erysipelotrichi";"Erysipelotrichales";Erysipelotrichaceae;Erysipelotrichaceae Incertae Sedis
376	0	0	0	0	0	0	0	0	1	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales
377	0	0	0	0	0	0	0	1	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Ruminococcaceae"
378	0	0	0	0	0	0	0	0	1	Root;Bacteria;Bacteroidetes
379	0	0	0	0	0	1	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Ruminococcaceae"
380	0	0	0	0	0	0	0	0	1	Root;Bacteria;Firmicutes;"Bacilli";Bacillales;"Staphylococcaceae";Staphylococcus
381	0	0	2	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
382	0	0	0	0	0	0	0	1	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales
383	4	9	0	2	0	0	0	2	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
384	0	0	1	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales
385	0	0	0	0	0	0	0	0	1	Root;Bacteria;Firmicutes;"Bacilli";"Lactobacillales";"Carnobacteriaceae";"Carnobacteriaceae 1"
386	0	0	1	0	0	0	0	0	0	Root;Bacteria
387	0	0	1	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
388	0	0	0	1	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
389	0	1	0	0	0	0	0	0	0	Root;Bacteria
390	0	0	0	0	0	0	0	0	1	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales
391	0	0	0	0	0	0	0	0	1	Root;Bacteria;Firmicutes
392	0	1	0	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales
393	0	0	0	0	0	1	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales
394	0	0	1	0	0	0	0	0	0	Root;Bacteria
395	1	1	1	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales
396	2	0	0	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
397	0	0	0	0	0	0	0	1	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
398	0	0	0	0	0	0	0	1	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
399	0	0	0	0	0	0	13	0	0	Root;Bacteria;Bacteroidetes;Bacteroidetes;Bacteroidales;Bacteroidaceae;Bacteroides
400	0	0	0	0	0	0	1	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
401	0	1	0	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
402	0	1	0	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
403	0	0	0	0	0	0	0	1	0	Root;Bacteria;Bacteroidetes;Bacteroidetes;Bacteroidales;Prevotellaceae
404	0	0	0	0	0	0	0	1	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae";"Lachnospiraceae Incertae Sedis"
405	0	0	0	0	0	0	0	1	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
406	0	0	0	0	0	1	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales
407	1	0	0	0	0	4	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales
408	1	5	3	2	0	0	0	0	1	Root;Bacteria;Bacteroidetes
409	0	0	0	0	0	0	0	1	1	Root;Bacteria;Bacteroidetes
410	0	0	0	0	1	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
411	0	0	0	1	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
412	0	0	0	0	2	0	0	0	0	Root;Bacteria;Bacteroidetes
413	0	0	0	0	0	0	0	1	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales
414	1	0	1	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
415	0	0	0	0	0	7	0	2	2	Root;Bacteria;Bacteroidetes
416	0	1	0	0	1	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales"""


if __name__ == "__main__":
    main()