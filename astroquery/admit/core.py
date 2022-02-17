# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""Download ADMIT data"""
import re
import os
# import path
import sqlite3
from astropy.table import Table
from ..query import BaseQuery
from ..utils import commons, async_to_sync

__all__ = ['ADMIT', 'ADMITClass']


@async_to_sync
class ADMITClass(BaseQuery):
    """
    TODO: document
    """

    request_url = 'http://admit.astro.umd.edu/query/'
    timeout = 60

    q = None
    db = None
    c  = None
    #  uid|ra|dec|z|object|
    rdz = 'datasetid_ra_dec_redshift_resolvername.txt'
    rdz_lines = None

    def query(self, *args, **kwargs):
        """
        query ADMIT
        """
        print("args   ",args)
        print("kwargs ",kwargs)
        if 'ADMIT' in os.environ:
            self.q = os.environ['ADMIT'] + '/query'
            self.db = self.q + '/admit.db'
            if os.path.exists(self.db):
                print("Found ",self.db)
                self.c = sqlite3.connect(self.db)
                print('Checking db....',self.c.total_changes)
            else:
                print("Did not find ",self.db)
        else:
            print("$ADMIT not in environment. Expecting $ADMIT/query/admit.db")

    def check(self):
        """
        """
        if self.c == None:
            print("database not open yet")
        for t in ["spw", "lines", "sources", "cont"]:
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
