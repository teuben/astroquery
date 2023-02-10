# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""Download ADMIT data"""
import re
import os
import sqlite3
import pickle
from astropy.table import Table
import pandas as pd
from astropy import units as u
from astropy.coordinates import SkyCoord
#from astropy.time import Time

from ..query import BaseQuery
from ..utils import commons, async_to_sync
from ..alma.tapsql import _gen_pos_sql, _gen_str_sql, _gen_numeric_sql,\
    _gen_band_list_sql, _gen_datetime_sql, _gen_pol_sql, _gen_pub_sql,\
    _gen_science_sql, _gen_spec_res_sql, ALMA_DATE_FORMAT

__all__ = ['ADMIT', 'ADMITClass','ADMIT_FORM_KEYS']
        

__version__ = "26-feb-2022"

def version():
    return __version__

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
        'Spectral window ID': [ 'win_id','win.id', _gen_numeric_sql],
        'Spectral window name': [ 'spw','win.spw', _gen_str_sql],
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
        'Spectral window ID(L)': [ 'l_win_id','lines.w_id', _gen_numeric_sql],# user should never select on this
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
        'Spectral Window ID(S)': [ 's_win_id','sources.w_id', _gen_numeric_sql], # user should never select on this
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
        'Search Region':['region','sources.region',None],# special case where we will parse internally
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
        # Energyi
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
        'Project abstract': ['project_abstract', 'alma.project_abstract', _gen_str_sql],
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
    a  = None     # alma pandas
    c  = None     # sqlite3 conn to admit.db
    ca = None     # sqlite3 conn to alma.db
    
    #  uid|ra|dec|z|object|
    rdz = 'datasetid_ra_dec_redshift_resolvername.txt'
    rdz_lines = None

    def __init__(self,db=None,pickle=False, alma=False,check_same_thread=True):
        '''if $ADMIT then we will look for $ADMIT/query/admit.db, otherwise db is full path to database to read '''
        # use check_same_thread=False if in a web app or you get an error.
        self.check_same_thread = check_same_thread
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
        print(self.report_version())
    
    def report_version(self):
        return f"Database version: {self.database_version}. core.py version: {version()}"

    @property
    def database_version(self):
        h = self.sql("select * from header;")
        hd = dict()
        for kv in h:
            hd[kv[1]] = kv[2]
        return hd['version']
    
    def load_admit(self, admit_db):
        '''Load the local ALMA+ADMIT database.  This method is necessary because the on-line ALMA archive does
        not yet include ADMIT products.
        
        Parameters:
            admit_db - str. Fully qualified path to database file
        '''
        if os.path.exists(admit_db):
            print("Found ",admit_db)
            self.c = sqlite3.connect(admit_db,
                                     check_same_thread=self.check_same_thread)
            print('Checking db....',self.c.total_changes)
        else:
            print("Did not find ",admit_db)

    def load_alma(self, alma_db):
        '''Load the local pure ALMA database.  This is a fallback method, only for debugging
        since the ADMIT database include an (albeit smaller) ALMA table.
        
        Parameters:
            alma_db - str. Fully qualified path to database file
        '''
        if os.path.exists(alma_db):
            print("Found ",alma_db)
            self.ca = sqlite3.connect(alma_db)
            print('Checking db....',self.ca.total_changes)
            print("Found %d entries" % len(self.sqla("SELECT obs_id from alma")))
            
        else:
            print("Did not find ",alma_db)
            
    def _set_colnames(self):
        '''Build the keyword dict and descriptive table'''
        type_dict = {'integer':int, 'float':float, 'text':str} #convert SQL to python type
        self._colnames = dict()
        self._coltypes = dict()
        for tab in ['alma','win','sources','lines','header']:
            result = self.sql(f"PRAGMA table_info({tab});")
            self._colnames[tab] = [x[1] for x in result]
            # _coltypes not currently used but perhaps helpful
            self._coltypes[tab] = [type_dict[x[2].lower()] for x in result]
            
    def _set_keys(self):
        #probably a more pythonic way to do this
        self.keys = list()
        self._key_description = dict()
        dbtable = dict()
        for attrib_category_key in ADMIT_FORM_KEYS:
            attrib_category = ADMIT_FORM_KEYS[attrib_category_key]
            for attrib_k in attrib_category:
                #self.keys.append(attrib[1].split(".")[1])
                self.keys.append(attrib_category[attrib_k][0])
                self._key_description[attrib_k] = attrib_category[attrib_k][0]
                dbtable[attrib_k] = attrib_category_key
        c1 = list(self._key_description.values())
        c2 = list(self._key_description.keys())
        c3 = list(dbtable.values())
        print(len(c1),len(c2),len(c3))
        self.ktable = Table([c1,c2,c3],names=("Keyword", "Description","Database Table"))
     
    def load_alma_pickle(self, alma_pickle):
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
        '''A table listing the possible keywords for query(), their descriptions, and what subtable they are in.'''
        return self.ktable
    
    def query(self, **kwargs):
        """
        Submit a query to the ALMA+ADMIT database. 
        Parameters:
            kwargs - string values for any keyword available in ADMIT.keys.
            Values MUST be strings, even if the underlying table column is number.  e.g.
            mom0flux='>0'  not mom0flux>0
        Returns: pandas DataFrame containing the search results table.
        """
        #print("kwargs ",kwargs)
        # alma and win are always present
        self._out_colnames = self._colnames['alma'] + self._colnames['win'] 
        self._out_coltypes = self._coltypes['alma'] + self._coltypes['win'] 
        if len(kwargs) == 0:
            raise Exception("You must supply at least one search keyword")
        bad = set(list(kwargs.keys())) - set(self.keys)
        if(len(bad)>0):
            raise Exception(f"Unrecognized keywords: {bad}")
        else:
            sqlq = self._gen_sql(kwargs)
            print(sqlq,"\n")
            return pd.DataFrame(data=self.sql(sqlq),columns=self._out_colnames)

    def check(self):
        """
        Internal check
        """
        if self.c == None:
            print("database not open yet")
        for t in ["header", "alma", "win", "lines", "sources"]:
            print("%-10s: %d entries" % (t,len(self.sql("SELECT id from %s" % t))))

    def sql(self, command):
        """
        Execute an arbitrary SQL; use with caution.
        Only meant for debugging
        Returns: list of lists of SQL result.
        """
        if False:
            cur = self.c.cursor()
            rows = cur.execute(command).fetchall()
        else:
            rows = self.c.execute(command).fetchall()
        return rows

    def sqla(self, command):
        """
        Execute an arbitrary SQL on the pure ALMA table; use with caution.
        Only meant for debugging
        For other experiments the alma.pickle approach, which contains
        the alma table as a panda's DataFrame, may be more efficient.
        Plus, there is always
            pd.read_sql('SELECT * FROM alma', self.ac)
        Returns: list of lists of SQL result.
        """
        return self.ca.execute(command).fetchall()

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
        
    def _parse_region(self,region):
        ''' get a square region centered on ra, dec '''
        c = region[0]
        size = region[1]
        if not isinstance(c,SkyCoord):
            raise ValueError("RA_DEC region value must be astropy SkyCoord")
        if not isinstance(size,u.Quantity):
            raise ValueError("SIZE region value must be astropy Quantities")
        ra_min = c.ra - size.to("degree")
        ra_max = c.ra + size.to("degree")
        dec_min = c.dec - size.to("degree")
        dec_max = c.dec + size.to("degree")
        sql = f" sources.ra > {ra_min.to('degree').value} AND sources.ra < {ra_max.to('degree').value} AND sources.dec > {dec_min.to('degree').value} AND sources.dec < {dec_max.to('degree').value} "
        return sql  
    
    def _gen_sql(self,payload):
        '''Transform the user keywords into an SQL string'''
        # Query joins alma and win tables always.
        sql = 'select * from alma inner join win on (win.a_id = alma.id) '
        # We may join sources table or lines table depending on what keywords the user invokes
        needs_lines_join = False
        needs_source_join = False 
        join_s = ' inner join sources on (sources.w_id = win.id)  '
        join_l = ' inner join lines on (lines.w_id = win.id ) '
        # Default is no additional joins
        join =  ''
        where = ''
        if payload:
            for constraint in payload:
                for attrib_category in ADMIT_FORM_KEYS.values():
                    for attrib in attrib_category.values():
                        if constraint in attrib:
                            # Set triggers for additional joins as needed
                            if attrib in ADMIT_FORM_KEYS["Sources"].values():
                                needs_source_join = True
                            if attrib in ADMIT_FORM_KEYS["Lines"].values():
                                needs_lines_join = True                     
                            # use the value and the second entry in attrib which
                            # is the new name of the column
                            val = payload[constraint]
                            #print(constraint)
                            # em_resolution is special-cased. see astroquery/alma/core.py
                            # i think this is a bug in alma/core.py should be 'spectral_resolution'
                            if constraint == 'alma.em_resolution': 
                                # em_resolution does not require any transformation
                                attrib_where = _gen_numeric_sql(constraint, val)
                            elif constraint == "region":
                                attrib_where = self._parse_region(val)
                            else:
                                attrib_where = attrib[2](attrib[1], val)
                            # Example of ADMIT virtual keyword
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
        ############################################################################
        # Do whatever additional joins are ncessary based on the triggers above.
        # The order of the join matters so that we get the column names correct and
        # the user always gets them in the same order.
        # We should always join in this order:
        # ALMA, WIN, [SOURCES], [LINES]
        ############################################################################
        if needs_source_join:
            join += join_s 
            self._out_colnames += self._colnames['sources']
            self._out_coltypes = self._coltypes['sources']
        if needs_lines_join:
            join += join_l 
            self._out_colnames += self._colnames['lines']
            self._out_coltypes = self._coltypes['lines']
            if needs_source_join: 
                # Do not move this to the if needs_sources_join above.
                # Order matters.  We should never have both sources.l_id=0 and sources.l_id>0
                where += ' AND sources.l_id > 0 ' #include cubesum results
        else:
            # Do not move this to the if needs_sources_join above. 
            # Order matters.  We should never have both sources.l_id=0 and sources.l_id>0
            if needs_source_join:
                where += ' AND sources.l_id = 0 ' #include linecube results
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


