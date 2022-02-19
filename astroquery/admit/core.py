# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""Download ADMIT data"""
import re
import os
import sqlite3
import pickle
from astropy.table import Table
from ..query import BaseQuery
from ..utils import commons, async_to_sync
from ..alma.tapsql import _gen_pos_sql, _gen_str_sql, _gen_numeric_sql,\
    _gen_band_list_sql, _gen_datetime_sql, _gen_pol_sql, _gen_pub_sql,\
    _gen_science_sql, _gen_spec_res_sql, ALMA_DATE_FORMAT

__all__ = ['ADMIT', 'ADMITClass']

# This mimics the ALMA_FORM_KEYS in alma/core.py.  The assumption here is
# that there is a web form in front of this.  We don't have it for ADMIT, but
# pretending we do means we can follow alma/core.py functions.
# format is 
# 'Category' : {
#  'Attribute1' :[ payload_kw, table_colname, function-for-parsing]
#  'Attribute2' :[ payload_kw, table_colname, function-for-parsing]
# }, 
#
ADMIT_FORM_KEYS = {
    'Window': {  
    # we need "spw.id and sources.spw_id" so put spw in two
        'ALMA ID': ['alma_id','spw.alma_id',_gen_numeric_sql],
        'Spectral Window': [ 'spw','spw.id', _gen_numeric_sql],
        'Number of Lines': [ 'nlines','spw.nlines', _gen_numeric_sql],
        'Number of Sources': [ 'nsources','spw.nsources', _gen_numeric_sql],
        'Number of Channels': [ 'nchan','spw.nchan', _gen_numeric_sql],
        'RMS noise': [ 'rms','spw.rms', _gen_numeric_sql],
     },
    'Lines': {
        'Spectral Window': [ 'spw','lines.spw_id', _gen_numeric_sql],
        'Transition': [ 'transition','lines.transition', _gen_str_sql],
        'Velocity': [ 'velocity','lines.velocity', _gen_numeric_sql],
        # we are not using this
        #'Channels': [ 'chan','lines.chan', _gen_numeric_sql],
     },
    'Sources': {
        'Spectral Window': [ 'spw','sources.spw_id', _gen_numeric_sql],
        'Line ID': [ 'lines_id','sources.lines_id', _gen_numeric_sql],
        'RA (Degrees)': ['ra', 'sources.ra',  _gen_numeric_sql],
        'Dec (Degrees)': ['dec', 'sources.dec',  _gen_numeric_sql],
        'Flux': ['flux', 'sources.flux',  _gen_numeric_sql],
     },
    'Header': { # no science use case
        'Key': ['header_key','header.key',_gen_str_sql],
        'Value': ['header_val','header.val',_gen_str_sql],
     },
    'Cont': {
        'ALMA ID': ['alma_id','cont.alma_id',_gen_numeric_sql],
        'Bands': ['band','cont.cont',_gen_str_sql],
        'Number of Continuum Sources': ['nc_sources','cont.nsources',_gen_numeric_sql],
     },
    'Alma': {
        'Observation': ['obs_id','alma.obs_id',_gen_str_sql],
        'Source name': ['source_name','alma.target_name',_gen_str_sql],
        'Source name': ['source_name','alma.target_name',_gen_str_sql],
        'RA Dec (Sexagesimal)': ['ra_dec', 'alma.s_ra, alma.s_dec', _gen_pos_sql],
        'Frequency (GHz)': ['frequency', 'alma.frequency', _gen_numeric_sql],
     },
}

@async_to_sync
class ADMITClass(BaseQuery):
    """
    TODO: document
    """

    request_url = 'http://admit.astro.umd.edu/query/'
    timeout = 60

    q  = None     # root query dir
    db = None     # admit db name
    pa = None     # alma pickle name
    c  = None     # sqlite3 conn
    a  = None     # alma pandas
    
    #  uid|ra|dec|z|object|
    rdz = 'datasetid_ra_dec_redshift_resolvername.txt'
    rdz_lines = None

    def __init__(self):
        if 'ADMIT' in os.environ:
            self.q = os.environ['ADMIT'] + '/query'
            # admit sqlite3
            self.db = self.q + '/admit.db'
            self.load_admit(self.db)
            # alma pickle
            self.pa = self.q + '/alma.pickle'
            self.load_alma(self.pa)
        else:
            print("$ADMIT not in environment. Expecting $ADMIT/query/admit.db")

    def load_admit(self, admit_db):
        if os.path.exists(admit_db):
            print("Found ",admit_db)
            self.c = sqlite3.connect(admit_db)
            print('Checking db....',self.c.total_changes)
        else:
            print("Did not find ",admit_db)
        
    def load_alma(self, alma_pickle):
        if os.path.exists(alma_pickle):
            print("Found ",alma_pickle)
            self.a = pickle.load(open(alma_pickle,'rb'))
            print("ALMA: Found %d entries" % len(self.a))
        else:
            print("Did not find ",alma_pickle)

    def query(self, *args, **kwargs):
        """
        query ADMIT
        """
        print("args   ",args)
        print("kwargs ",kwargs)
        print(_gen_sql(payload))

    def check(self):
        """
        """
        if self.c == None:
            print("database not open yet")
        for t in ["header", "alma", "spw", "lines", "sources", "cont"]:
            print("%-10s: %d entries" % (t,len(self.sql("SELECT id from %s" % t))))

    def sql(self, command):
        """
        execute an arbitrary SQL; use with caution
        only meant for debugging
        """
        if False:
            cur = self.c.cursor()
            rows = cur.execute(command).fetchall()
        else:
            rows = self.c.execute(command).fetchall()
        return rows

    def query_rdz(self, *args, **kwargs):
        """
        query the ALMA ra-dec-redshift catalog (courtesy: F. Stoehr)
        """
        if self.rdz_lines == None:
            self.rdz_lines = open(self.q + '/datasetid_ra_dec_redshift_resolvername.txt').readlines()
            print("Found %d entries in rdz" % len(self.rdz_lines))

    def query_sql_async(self, *args, **kwargs):
        """
        Query the ADMIT database

        Returns
        -------
        url : The URL of the FITS file containing the results.
        """

        payload = self._parse_args(*args, **kwargs)

        if kwargs.get('get_query_payload'):
            return payload

        result = self._request("POST", url=self.request_url,
                               data=payload, timeout=self.timeout)

        result_url_relative = find_data_url(result.text)
        result_url = os.path.join(self.request_url, result_url_relative)

        return result_url

    def _parse_args(self, sql_query):
        """
        Parameters
        ----------
        sql_query : str
            An SQL query

        Returns
        -------
        payload_dict : Requests payload in a dictionary
        """

        payload = {'query': sql_query,
                   'format': 'fits'}

        return payload

    def _parse_result(self, result, verbose=False, **kwargs):
        """
        Use get_admit_datafile to download a result URL
        """
        return get_admit_datafile(result)

    def help(self):
        """
        reminders...
        """
        print("admit.db contains tables that work with the ALMA table")
        print("   alma:      a placeholder for the big alma table")
        print("   spw:       all SPW's are stored here, with links back to the alma reference")
        print("   lines:     all lines detected in the SPW's")
        print("   sources:   all sources detected in each ADMIT LineCube")


ADMIT = ADMITClass()


def get_admit_datafile(result, **kwargs):
    """Turn a URL into an HDUList object."""
    fitsfile = commons.FileContainer(result,
                                     encoding='binary',
                                     **kwargs)
    hdulist = fitsfile.get_fits()
    return Table(hdulist[1].data)


def find_data_url(result_page):
    """Find and return the URL of the data, given a results page."""
    result_relative_url_re = re.compile(r'Download the result file: '
                                        r'<a href="(\.\./tmp/.*?)">')
    re_result = result_relative_url_re.findall(result_page)
    if len(re_result) == 0:
        raise ValueError("Results did not contain a result url")
    return re_result[0]


# maybe we will have a different sql param for inner join
def _gen_sql2(payload,sql = 'select * spw,lines,sources'):
    pass

def _gen_sql(payload):
    sql = 'select * spw,lines,sources'
    where = ''
    if payload:
        for constraint in payload:
            for attrib_category in ADMIT_FORM_KEYS.values():
                for attrib in attrib_category.values():
                    if constraint in attrib:
                        # use the value and the second entry in attrib which
                        # is the new name of the column
                        val = payload[constraint]
                        attrib_where = attrib[2](attrib[1], val)
                        if attrib_where:
                            if where:
                                where += ' AND '
                            else:
                                where = ' WHERE '
                            where += attrib_where
    return sql + where
