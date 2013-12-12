#!/usr/bin/env python
# File created on 16 Feb 2011
from __future__ import division

__author__ = "Jesse Stombaugh"
__copyright__ = "Copyright 2011, The QIIME Project"
__credits__ = ["Jesse Stombaugh"]
__license__ = "GPL"
__version__ = "1.8.0-dev"
__maintainer__ = "Jesse Stombaugh"
__email__ = "jesse.stombaugh@colorado.edu"

from urllib2 import urlopen, URLError
from qiime.util import split_fasta_on_sample_ids_to_files
from cogent.parse.fasta import MinimalFastaParser
import re
import os
from glob import glob
import urllib
import httplib
import mimetypes

log_html = '''\
<html>\
<head><title>QIIME to MG-RAST submission</title></head>\
<body>
%s\
</body>\
</html>\
'''


def parse_and_submit_params(key, project_id, seq_file, output_dir,
                            submit_to_server=True):
    '''This function takes the input options from the user and generates a url
       and request header for submitting to the MG-RAST cgi script'''

    # Verify that the users computer can connect to the internet
    try:
        check_internet = urlopen('http://www.google.com')
    except:
        raise OSError(
            "This script is having trouble connecting to the internet!")

    # parse and split fasta file into individual sample fastas
    fasta_file = MinimalFastaParser(open(seq_file))
    split_fasta_on_sample_ids_to_files(fasta_file, output_dir)

    # set the MG-RAST link for QIIME
    host = 'metagenomics.anl.gov'

    # open the log html
    log_file = open(os.path.join(output_dir, 'log.html'), 'w')
    log_data = ['<h3>The following jobs were submitted to MG-RAST.</h3>']
    log_data.append('<table border=1><tr><th>Fasta File</th><th>Job ID</th>')
    log_data.append('<th>md5</th></tr>')
    num = 0
    # iterate over the fasta files in the given directory
    fasta_filepaths = sorted(glob('%s/*.fasta' % output_dir))
    for i in fasta_filepaths:

        # Get the sample id from the fasta filename
        sample_id = os.path.split(os.path.splitext(i)[0])[-1]

        # set the parameters
        params = [('key', key), ('sample', sample_id), ('project', project_id)]

        # get the full path and short name for the fasta file to be uploaded
        file_to_submit = os.path.abspath(i)
        fasta_shortname = os.path.split(file_to_submit)[-1]

        # open and read file to be put in post form
        file_object = open(file_to_submit).read()

        # set the file
        files = [('file', fasta_shortname, file_object)]

        # Post the file and parameters
        response = post_multipart(host, params, files, submit_to_server)

        # check the response for MG-RAST errors
        job = re.findall(r'<id>.*</id>', response)
        md5 = re.findall(r'<md5>.*</md5>', response)

        # if job successful write to log html otherwise post an error message
        # in the log file
        if job and md5:
            job_id = job[0].strip('<id>').strip('</id>')
            md5_id = md5[0].strip('<md5>').strip('</md5>')
            log_data.append('<tr><td>%s</td><td>%s</td><td>%s</td></tr>' %
                            (fasta_shortname, job_id, md5_id))
        else:
            response_error = re.findall(
                r'Can\'t call method "login" ',
                response)
            if response_error:
                log_data.append('</table><br><h3 style="color:red">')
                log_data.append('Web-service authorization key is not valid!')
                log_data.append('</h3>')
            else:
                log_data.append('</table><br><h3 style="color:red">%s</h3>' %
                                (response))

    log_data.append('</table>')

    log_info = '\n'.join(log_data)
    # write and close the log html
    log_file.write(log_html % (log_info))
    log_file.close()

    return log_info

# {{{ http://code.activestate.com/recipes/146306/ (r1)
"""This function was taken from a recipe on the activestate website
   http://code.activestate.com/recipes/146306/ where the function was adapted
   to the MG-RAST cgi script"""


def post_multipart(host, fields, files, submit_to_server):
    """
    Post fields and files to an http host as multipart/form-data.
    fields is a sequence of (name, value) elements for regular form fields.
    files is a sequence of (name, filename, value) elements for data to be
    uploaded as files
    Return the server's response page.
    """
    content_type, body = encode_multipart_formdata(fields, files)
    h = httplib.HTTP(host)
    # needed to change the following url to be handled properly by MG-RAST
    h.putrequest('POST', 'http://metagenomics.anl.gov/qiime.cgi')
    h.putheader('Content-Type', content_type)
    h.putheader('Content-Length', str(len(body)))
    h.endheaders()

    # put a check in place for testing purposes on whether the data should be
    # posted on the MG-RAST website
    if submit_to_server:
        h.send(body)
        errcode, errmsg, headers = h.getreply()

        # verify the data was received by MG-RAST
        if errcode == 200:
            response = h.file.read()
        else:
            raise OSError(
                'MG-RAST could not fulfill the request, which means that the server is unavailable!')
    else:
        response = body

    return response

"""This function was taken from a recipe on the activestate website
   http://code.activestate.com/recipes/146306/ where the function was adapted
   to the MG-RAST cgi script"""


def encode_multipart_formdata(fields, files):
    """
    fields is a sequence of (name, value) elements for regular form fields.
    files is a sequence of (name, filename, value) elements for data to be
    uploaded as files
    Return (content_type, body) ready for httplib.HTTP instance
    """
    # changed the boundary to be more similar to the perl script written by
    # Andreas
    BOUNDARY = 'xYzZY'
    CRLF = '\r\n'
    L = []
    for (key, value) in fields:
        L.append('--' + BOUNDARY)
        L.append('Content-Disposition: form-data; name="%s"' % key)
        L.append('')
        L.append(value)
    for (key, filename, value) in files:
        L.append('--' + BOUNDARY)
        L.append('Content-Disposition: form-data; name="%s"; filename="%s"' %
                 (key, filename))
        L.append('Content-Type: %s' % get_content_type(filename))
        L.append('')
        L.append(value)
    L.append('--' + BOUNDARY + '--')
    L.append('')
    body = CRLF.join(L)
    content_type = 'multipart/form-data'
    return content_type, body

"""This function was taken from a recipe on the activestate website
   http://code.activestate.com/recipes/146306/ where the function was adapted
   to the MG-RAST cgi script"""


def get_content_type(filename):
    # changed the default to text/plain since it is the format of fasta files
    return mimetypes.guess_type(filename)[0] or 'text/plain'

# end of http://code.activestate.com/recipes/146306/ }}}
