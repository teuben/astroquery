Contributing to Astroquery
==========================
The first official release of astroquery occurred in September 2013.

Please see `astropy's contributing guildelines
<https://www.astropy.org/contribute.html>`__ for a general guide to the
workflow involving git, etc.  Everything below is astroquery-specific.

We strongly encourage draft pull requests to be opened early in development.
If you are thinking of contributing a new module, please open a pull request
as soon as you start developing code and mark it as a Draft PR on github.


New Features
------------
We welcome any and all new features!  If you have your own little query tool
you wrote to access some obscure service, feel free to clean it up a little and
submit it as a pull request (PR)!  However, we'll ask you to at least write
some simple tests and do your best to get the API to match the `astroquery API`_.

Tests are welcome and encouraged, but can be built up over time.  At least one
example use is necessary for a new module to be accepted.

The template_ module should be the starting point for any new contribution.
It includes detailed explanations of each function and class and can
straightforwardly be filled out.

As shown in the template_ module, the minimum requirements for a new feature are:

 * Add the feature as a subdirectory of astroquery with at least an
   ``__init__.py`` and a ``core.py``::
 
     astroquery/feature
     astroquery/feature/__init__.py
     astroquery/feature/core.py

 * Add a ``tests/`` directory with at least one test::
 
     astroquery/feature/tests
     astroquery/feature/tests/__init__.py
     astroquery/feature/tests/test_feature.py

 * Add some documentation - at least one example, but it can be sparse at first::
 
     docs/astroquery/feature.rst

PRs will not be accepted if tests fail.

Important Guidelines
--------------------

Astroquery is based on the requests module.  All contributions must be based on
the `requests`_ module and should not use `urllib` or any of the base python url
modules unless there is a demonstrated necessity.

The `requests`_ module also generally should not be directly used, since the
`astroquery.query.BaseQuery` class, which all astroquery classes should inherit
from, provides access to its own `_request` method.  This custom `_request`
method is a wrapper around the `requests.request` function that provides
important astroquery-specific utility, including caching, HTTP header
generation, progressbars, and local writing-to-disk.

Dependencies
------------
New contributions are generally not allowed to bring along additional dependencies.

The astropy ecosystem tools should be used whenever possible.
For example, `astropy.table` should be used for table handling,
or `astropy.units` for unit and quantity
handling.



.. _astroquery API: docs/api.rst
.. _template: docs/template.rst
.. _requests: https://docs.python-requests.org/en/latest/
