from colabdesign.af.model import mk_af_model
from colabdesign.tr.model import mk_tr_model
from colabdesign.tr.joint_model import mk_af_tr_model
from colabdesign.shared.utils import clear_mem

# backward compatability
mk_design_model = mk_afdesign_model = mk_af_model
mk_trdesign_model = mk_tr_model