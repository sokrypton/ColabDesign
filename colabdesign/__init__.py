import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

from colabdesign.shared.utils import clear_mem
from colabdesign.af.model import mk_af_model
from colabdesign.tr.model import mk_tr_model
from colabdesign.mpnn.model import mk_mpnn_model

# backward compatability
mk_design_model = mk_afdesign_model = mk_af_model
mk_trdesign_model = mk_tr_model