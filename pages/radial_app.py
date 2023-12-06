'''
                    Waveplot: An online wavefunction viewer
                    Copyright (C) 2023  Jon G. C. Kragskow

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU Affero General Public License as published
    by the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Affero General Public License for more details.

    You should have received a copy of the GNU Affero General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.
'''

from dash import register_page
from .core import radial as rc
from .core import utils
from .core import common

ID_PREFIX = 'rad'
dash_id = common.dash_id(ID_PREFIX)

PAGE_NAME = 'Radial Wavefunctions'
PAGE_PATH = '/radial'
PAGE_IMAGE = 'assets/radial.png'
PAGE_DESCRIPTION = 'Interactive Radial Functions'

register_page(
    __name__,
    order=1,
    path=PAGE_PATH,
    name=PAGE_NAME,
    title=PAGE_NAME,
    image=PAGE_IMAGE,
    description=PAGE_DESCRIPTION
)

prefix = utils.dash_id('orb3')

plot_div = rc.PlotDiv(ID_PREFIX)

# Make AC options tab and all callbacks
options = rc.OptionsDiv(ID_PREFIX)
# Connect callbacks for plots and options
rc.assemble_callbacks(plot_div, options)

# Layout of webpage
layout = common.make_layout(plot_div.div, options.div)