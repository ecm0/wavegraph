# -*- coding: utf-8 -*-
# Copyright (C) Eric Chassande-Mottin (2014-2019)
#
# This file is part of pyburst
#
# pyburst is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# pyburst is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with pyburst.  If not, see <http://www.gnu.org/licenses/>.

"""A Python package for searching gravitational-wave burst signals
"""

import logging
LOGGING_FMT = "%(levelname)s -- %(filename)s:line %(lineno)s in %(funcName)s(): %(message)s"

from ._version import get_versions
__version__ = get_versions()['version']
__author__ = "Eric Chassande-Mottin <eric.chassandemottin@ligo.org>"
__credits__ = "The LIGO Scientific Collaboration and the Virgo Collaboration"

del get_versions
