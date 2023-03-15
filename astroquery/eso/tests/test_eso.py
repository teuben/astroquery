# Licensed under a 3-clause BSD style license - see LICENSE.rst
import os

from astroquery.utils.mocks import MockResponse
from ...eso import Eso

DATA_DIR = os.path.join(os.path.dirname(__file__), 'data')


def data_path(filename):
    return os.path.join(DATA_DIR, filename)


DATA_FILES = {
    'GET':
        {
            'http://archive.eso.org/wdb/wdb/eso/eso_archive_main/form': 'main_query_form.html',
            'http://archive.eso.org/wdb/wdb/eso/amber/form': 'amber_query_form.html',
            'http://archive.eso.org/wdb/wdb/adp/phase3_main/form': 'vvv_sgra_form.html',
            Eso.AUTH_URL: 'oidc_token.json',
        },
    'POST':
        {
            'http://archive.eso.org/wdb/wdb/eso/eso_archive_main/query': 'main_sgra_query.tbl',
            'http://archive.eso.org/wdb/wdb/eso/amber/query': 'amber_sgra_query.tbl',
            'http://archive.eso.org/wdb/wdb/adp/phase3_main/query': 'vvv_sgra_survey_response.tbl',
        }
}


def eso_request(request_type, url, **kwargs):
    with open(data_path(DATA_FILES[request_type][url]), 'rb') as f:
        response = MockResponse(content=f.read(), url=url)
    return response


def download_request(url, **kwargs):
    filename = 'testfile.fits.Z'
    with open(data_path(filename), 'rb') as f:
        header = {'Content-Disposition': f'filename={filename}'}
        response = MockResponse(content=f.read(), url=url, headers=header)
    return response


def calselector_request(url, **kwargs):
    is_multipart = len(kwargs['data']['dp_id']) > 1
    if is_multipart:
        filename = 'FORS2.2021-01-02T00:59:12.533_raw2raw_multipart.xml'
        header = {
            'Content-Type': 'multipart/form-data; boundary=uFQlfs9nBIDEAIoz0_ZM-O2SXKsZ2iSd4h7H;charset=UTF-8'
        }
    else:
        filename = 'FORS2.2021-01-02T00:59:12.533_raw2raw.xml'
        header = {
            'Content-Disposition': f'filename="{filename}"',
            'Content-Type': 'application/xml; content=calselector'
        }
    with open(data_path(filename), 'rb') as f:
        response = MockResponse(content=f.read(), url=url, headers=header)
    return response


# @pytest.fixture
# def patch_get(request):
#    mp = request.getfixturevalue("monkeypatch")
#
#    mp.setattr(Eso, 'request', eso_request)
#    return mp

# This test should attempt to access the internet and therefore should fail
# (_activate_form always connects to the internet)
# @pytest.mark.xfail
def test_amber_SgrAstar(monkeypatch):
    # Local caching prevents a remote query here

    eso = Eso()

    # monkeypatch instructions from https://pytest.org/latest/monkeypatch.html
    monkeypatch.setattr(eso, '_request', eso_request)
    # set up local cache path to prevent remote query
    eso.cache_location = DATA_DIR

    # the failure should occur here
    result = eso.query_instrument('amber', target='Sgr A*')

    # test that max_results = 50
    assert len(result) == 50
    assert 'GC_IRS7' in result['Object']


def test_main_SgrAstar(monkeypatch):
    # Local caching prevents a remote query here

    eso = Eso()

    # monkeypatch instructions from https://pytest.org/latest/monkeypatch.html
    monkeypatch.setattr(eso, '_request', eso_request)
    # set up local cache path to prevent remote query
    eso.cache_location = DATA_DIR

    # the failure should occur here
    result = eso.query_main(target='Sgr A*')

    # test that max_results = 50
    assert len(result) == 50
    assert 'GC_IRS7' in result['OBJECT']


def test_vvv(monkeypatch):
    eso = Eso()
    monkeypatch.setattr(eso, '_request', eso_request)
    eso.cache_location = DATA_DIR

    result_s = eso.query_surveys(surveys='VVV',
                                 coord1=266.41681662, coord2=-29.00782497,
                                 box='01 00 00',
                                 )
    assert result_s is not None
    assert 'Object' in result_s.colnames
    assert 'b333' in result_s['Object']


def test_authenticate(monkeypatch):
    eso = Eso()
    monkeypatch.setattr(eso, '_request', eso_request)
    eso.cache_location = DATA_DIR
    authenticated = eso._authenticate(username="someuser", password="somepassword")
    assert authenticated is True


def test_download(monkeypatch):
    fileid = 'testfile'
    url = Eso.DOWNLOAD_URL + fileid
    eso = Eso()
    destination = os.path.join(DATA_DIR, 'downloads')
    os.makedirs(destination, exist_ok=True)
    monkeypatch.setattr(eso._session, 'get', download_request)
    filename, downloaded = eso._download_eso_file(url, destination=destination, overwrite=True)
    assert downloaded is True
    assert fileid in filename


def test_calselector(monkeypatch):
    eso = Eso()
    dataset = 'FORS2.2021-01-02T00:59:12.533'
    monkeypatch.setattr(eso._session, 'post', calselector_request)
    result = eso.find_associated_files([dataset], savexml=True, destination=data_path('downloads'))
    assert isinstance(result, list)
    assert len(result) == 50
    assert dataset not in result


def test_calselector_multipart(monkeypatch):
    eso = Eso()
    datasets = ['FORS2.2021-01-02T00:59:12.533', 'FORS2.2021-01-02T00:59:12.534']
    monkeypatch.setattr(eso._session, 'post', calselector_request)
    result = eso.find_associated_files(datasets, savexml=True, destination=data_path('downloads'))
    assert isinstance(result, list)
    assert len(result) == 99
    assert datasets[0] not in result and datasets[1] not in result
