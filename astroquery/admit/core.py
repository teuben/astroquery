# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""Download ADMIT data"""
import re
import os
import sqlite3
import pickle
from astropy.table import Table
import pandas as pd
#from astropy import units as u
#from astropy.time import Time

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
        'ALMA ID': ['alma_id','win.alma_id',_gen_numeric_sql],
        'Spectral window': [ 'spw','win.id', _gen_numeric_sql],
        'Number of lines': [ 'nlines','win.nlines', _gen_numeric_sql],
        'Number of sources': [ 'nsources','win.nsources', _gen_numeric_sql],
        'Number of channels': [ 'nchan','win.nchan', _gen_numeric_sql],
        'Window peak flux': [ 'win_peak','win.peak_w', _gen_numeric_sql],
        'Window RMS noise': [ 'win_rms','win.rms', _gen_numeric_sql],
        'Window Signal to noise ratio (peak/rms)': [ 'win_snr','win.snr_w', _gen_numeric_sql],
        'Beam major axis':['bmaj','win.bmaj',_gen_numeric_sql],    
        'Beam minor axis':['bmin','win.bmin',_gen_numeric_sql],    
        'Beam PA':['bpa','win.bpaj',_gen_numeric_sql],  
        'Frequency center (GHz)':['freqc','win.freqc',_gen_numeric_sql],   
        'Frequency width (GHz)':['freqw','win.freqw',_gen_numeric_sql],  
        'LSR Velocity (km/s)':['vlsr','win.vlsr',_gen_numeric_sql],
        'Frequency coverage?':['fcoverage','win.fcoverage',_gen_numeric_sql], 
     },
    'Lines': {
        'Spectral window': [ 'spw','lines.w_id', _gen_numeric_sql],
        'Rest frequency': [ 'restfreq','lines.restfreq', _gen_numeric_sql],
        'Formula': [ 'formula','lines.formula', _gen_str_sql],
        'Transition': [ 'transition','lines.transition', _gen_str_sql],
        'Velocity': [ 'velocity','lines.velocity', _gen_numeric_sql],
        # would be good to have FWHM or virtual linewidth keyword that is vmax-vmin or something
        'Minimum velocity': ['vmin','lines.vmin', _gen_numeric_sql],   
        'Maximum velocity': ['vmax','lines.vmin', _gen_numeric_sql],
        'Moment zero flux': ['mom0flux','lines.mom0flux', _gen_numeric_sql],#Jy km/s?
        'Moment one peak': ['mom1peak','lines.mom1peak', _gen_numeric_sql],
        'Moment two peak (km/s)': ['mom2peak','lines.mom2peak', _gen_numeric_sql], 
    },
    'Sources': {
        'Spectral Window': [ 'spw','sources.w_id', _gen_numeric_sql],
        'Line ID': [ 'lines_id','sources.lines_id', _gen_numeric_sql],
        'RA (Degrees)': ['ra', 'sources.ra',  _gen_numeric_sql],
        'Dec (Degrees)': ['dec', 'sources.dec',  _gen_numeric_sql],
        'Flux': ['flux', 'sources.flux',  _gen_numeric_sql],
        'Source Signal to noise ratio': ['source_snr', 'sources.snr_s',  _gen_numeric_sql],
        'Source Peak flux':['source_peak','sources.peak_s',_gen_numeric_sql],
        # this cries out for a source size virtual keyword that is geometric mean of smaj&smin
        'Source major axis':['smaj','sources.smaj',_gen_numeric_sql],    
        'Source minor axis':['smin','sources.smin',_gen_numeric_sql],    
        'Source PA':['spa','sources.spaj',_gen_numeric_sql],  
     },
    'Header': { # no science use case
        'Key': ['header_key','header.key',_gen_str_sql],
        'Value': ['header_val','header.val',_gen_str_sql],
     },
#    'Cont': {
#        'ALMA ID': ['alma_id','cont.alma_id',_gen_numeric_sql],
#        'Bands': ['band','cont.cont',_gen_str_sql],
#        'Number of continuum sources': ['nc_sources','cont.nsources',_gen_numeric_sql],
#     },
    'Alma': {
        'Observation': ['obs_id','alma.obs_id',_gen_str_sql],
        # From here below are just a copy of ALMA_FORM_KEYS without the external wrapper dict.
        # Position
        'Source name (astropy Resolver)': ['source_name_resolver',
                                           'SkyCoord.from_name', _gen_pos_sql],
        'Source name (ALMA)': ['source_name_alma', 'alma.target_name', _gen_str_sql],
        'RA Dec (Sexagesimal)': ['ra_dec', 'alma.s_ra, alma.s_dec', _gen_pos_sql],
        'Galactic (Degrees)': ['galactic', 'alma.gal_longitude, alma.gal_latitude',
                               _gen_pos_sql],
        'Angular resolution (arcsec)': ['spatial_resolution',
                                        'alma.spatial_resolution', _gen_numeric_sql],
        'Largest angular scale (arcsec)': ['spatial_scale_max',
                                           'alma.spatial_scale_max', _gen_numeric_sql],
        'Field of view (arcsec)': ['fov', 'alma.s_fov', _gen_numeric_sql],
        # Energy
        'Frequency (GHz)': ['frequency', 'alma.frequency', _gen_numeric_sql],
        'Bandwidth (Hz)': ['bandwidth', 'alma.bandwidth', _gen_numeric_sql],
        'Spectral resolution (KHz)': ['spectral_resolution',
                                      'alma.em_resolution', _gen_spec_res_sql],
        'Band': ['band_list', 'alma.band_list', _gen_band_list_sql],
        # Time
        'Observation date': ['start_date', 'alma.t_min', _gen_datetime_sql],
        'Integration time (s)': ['integration_time', 'alma.t_exptime',
                                 _gen_numeric_sql],
        # Polarization
        'Polarisation type (Single, Dual, Full)': ['polarisation_type',
                                                   'alma.pol_states', _gen_pol_sql],
        # Observation
        'Line sensitivity (10 km/s) (mJy/beam)': ['line_sensitivity',
                                                  'alma.sensitivity_10kms',
                                                  _gen_numeric_sql],
        'Continuum sensitivity (mJy/beam)': ['continuum_sensitivity',
                                             'alma.cont_sensitivity_bandwidth',
                                             _gen_numeric_sql],
        'Water vapour (mm)': ['water_vapour', 'alma.pvw', _gen_numeric_sql],
        # Project
        'Project code': ['project_code', 'alma.proposal_id', _gen_str_sql],
        'Project title': ['project_title', 'alma.obs_title', _gen_str_sql],
        'PI name': ['pi_name', 'alma.obs_creator_name', _gen_str_sql],
        'Proposal authors': ['proposal_authors', 'alma.proposal_authors', _gen_str_sql],
        'Project abstract': ['project_abstract', 'alma.proposal_abstract', _gen_str_sql],
        'Publication count': ['publication_count', 'alma.NA', _gen_str_sql],
        'Science keyword': ['science_keyword', 'alma.science_keyword', _gen_str_sql],
        'Scientific category': ['scientific_category', 'alma.scientific_category', _gen_str_sql],
       # Publication'
        'Bibcode': ['bibcode', 'alma.bib_reference', _gen_str_sql],
        'Title': ['pub_title', 'alma.pub_title', _gen_str_sql],
        'First author': ['first_author', 'alma.first_author', _gen_str_sql],
        'Authors': ['authors', 'alma.authors', _gen_str_sql],
        # this may need special handling? or person does search on pub_abstract="*YSO*" 
        'Abstract': ['pub_abstract', 'alma.pub_abstract', _gen_str_sql],
        'Year': ['publication_year', 'alma.pub_year', _gen_numeric_sql],
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

    def __init__(self,db=None,pickle=False):
        '''if $ADMIT then we will look for $ADMIT/query/admit.db, otherwise db is full path to database to read '''
        if 'ADMIT' in os.environ and db is None:
            self.q = os.environ['ADMIT'] + '/query'
            # admit sqlite3
            self.db = self.q + '/admit.db'
            self.load_admit(self.db)
            # alma pickle, deprecated
            if pickle:
                self.pa = self.q + '/alma.pickle'
                self.load_alma(self.pa)
        else:
            if db is None:
                print("No db set yet - use load_admit(...)")
            else:
                self.db = db
                self.load_admit(self.db)
        self._set_keys()
        self._set_colnames()

    def load_admit(self, admit_db):
        if os.path.exists(admit_db):
            print("Found ",admit_db)
            self.c = sqlite3.connect(admit_db)
            print('Checking db....',self.c.total_changes)
        else:
            print("Did not find ",admit_db)
            
    def _set_colnames(self):
        type_dict = {'integer':int, 'float':float, 'text':str} #convert SQL to python type
        self._colnames = dict()
        self._coltypes = dict()
        for tab in ['alma','win','sources','lines','header']:
            result = self.sql(f"PRAGMA table_info({tab});")
            self._colnames[tab] = [x[1] for x in result]
            self._coltypes[tab] = [type_dict[x[2].lower()] for x in result]
            
    def load_alma(self, alma_pickle):
        """
        deprecated
        """
        if os.path.exists(alma_pickle):
            print("Found ",alma_pickle)
            self.a = pickle.load(open(alma_pickle,'rb'))
            print("ALMA: Found %d entries" % len(self.a))
        else:
            print("Did not find ",alma_pickle)
   
    @property
    def key_description(self):
        return self.ktable
    
    def _set_keys(self):
        #probably a more pythonic way to do this
        self.keys = list()
        self._key_description = dict()
        for attrib_category in ADMIT_FORM_KEYS.values():
            for attrib_k in attrib_category:
                #self.keys.append(attrib[1].split(".")[1])
                self.keys.append(attrib_category[attrib_k][0])
                self._key_description[attrib_k] = attrib_category[attrib_k][0]
        c1 = list(self._key_description.values())
        c2 = list(self._key_description.keys())
        self.ktable = Table([c1,c2],names=("Keyword", "Description"))

    def query(self, **kwargs):
        """
        query ADMIT
        """
        print("kwargs ",kwargs)
        # alma and win are always present
        self._out_colnames = self._colnames['alma'] + self._colnames['win'] + self._colnames['sources'] 
        self._out_coltypes = self._coltypes['alma'] + self._coltypes['win'] + self._coltypes['sources']
        if len(kwargs) == 0:
            raise Exception("You must supply at least one search keyword")
        bad = set(list(kwargs.keys())) - set(self.keys)
        if(len(bad)>0):
            print("WARNING: These keywords are not valid:", bad)
        else:
            sqlq = self._gen_sql(kwargs)
            print(sqlq)
            return pd.DataFrame(data=self.sql(sqlq),columns=self._out_colnames)
            #return self.sql(sqlq)

    def check(self):
        """
        """
        if self.c == None:
            print("database not open yet")
        for t in ["header", "alma", "win", "lines", "sources"]:
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
        print("   win:       all SPW's are stored here, with links back to the alma reference")
        print("   lines:     all lines detected in the SPW's")
        print("   sources:   all sources detected in each ADMIT LineCube")
        
    def _gen_sql(self,payload):
        needs_lines_join = False
        needs_source_join = False # I think we always want this join.
        sql = 'select * from alma inner join win on (win.a_id = alma.id) inner join sources on (sources.w_id = win.id)  '
        join =  ''
        where = ''
        if payload:
            for constraint in payload:
                for attrib_category in ADMIT_FORM_KEYS.values():
                    for attrib in attrib_category.values():
                        if constraint in attrib:
                            if attrib in ADMIT_FORM_KEYS["Sources"].values():
                                needs_source_join = True
                            if attrib in ADMIT_FORM_KEYS["Lines"].values():
                                needs_lines_join = True
                                needs_source_join = True     # if we have lines we also have sources                      
                            # use the value and the second entry in attrib which
                            # is the new name of the column
                            val = payload[constraint]
                            # see astroquery/alma/core.py
                            if constraint == 'alma.em_resolution':
                                # em_resolution does not require any transformation
                                attrib_where = _gen_numeric_sql(constraint, val)
                            else:
                                attrib_where = attrib[2](attrib[1], val)
                            # ADMIT virtual keyword
                            #  Replace win.snr_w with '( win.peak_ / win.rms_w )'
                            #  to compute signal to noise
                            if 'win.snr_w' in attrib_where:
                                attrib_where = attrib_where.replace('win.snr_w','( win.peak_w / win.rms_w )')
                            if attrib_where:
                                if where:
                                    where += ' AND '
                                else:
                                    where = ' WHERE '
                                where += attrib_where  
    #    if needs_source_join:
    #        join += ' inner join sources on (sources.w_id = win.id) '
        if needs_lines_join:
            join += ' inner join lines on (lines.w_id = win.id ) '
            self._out_colnames += self._colnames['lines']
            self._out_coltypes = self._coltypes['lines']
        if join:
            sql = sql + join

        return sql + where

ADMIT = ADMITClass


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


