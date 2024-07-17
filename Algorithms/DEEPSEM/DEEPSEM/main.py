import argparse

from src.DeepSEM_cell_type_non_specific_GRN_model import non_celltype_GRN_model
from src.DeepSEM_cell_type_specific_GRN_model import celltype_GRN_model
from src.DeepSEM_cell_type_test_non_specific_GRN_model import test_non_celltype_GRN_model
from src.DeepSEM_cell_type_test_specific_GRN_model import celltype_GRN_model as test_celltype_GRN_model
from src.DeepSEM_embed_model import deepsem_embed
from src.DeepSEM_generation_model import deepsem_generation

parser = argparse.ArgumentParser()
parser.add_argument('--n_epochs', type=int, default=120, help='Number of Epochs for training DeepSEM')
parser.add_argument('--task', type=str, default='celltype_GRN',
                    help='Determine which task to run. Select from (non_celltype_GRN,celltype_GRN,embedding,simulation)')
parser.add_argument('--setting', type=str, default='default', help='Determine whether or not to use the default hyper-parameter')
parser.add_argument('--batch_size', type=int, default=64, help='The batch size used in the training process.')
parser.add_argument('--data_file', type=str, help='The input scRNA-seq gene expression file.')
parser.add_argument('--net_file', type=str, default='',
                    help='The ground truth of GRN. Only used in GRN inference task if available. ')
parser.add_argument('--alpha', type=float, default=100, help='The loss coefficient for L1 norm of W, which is same as \\alpha used in our paper.')
parser.add_argument('--beta', type=float, default=1, help='The loss coefficient for KL term (beta-VAE), which is same as \\beta used in our paper.')
parser.add_argument('--lr', type=float, default=1e-4, help='The learning rate of used for RMSprop.')
parser.add_argument('--lr_step_size', type=int, default=0.99, help='The step size of learning rate decay.')
parser.add_argument('--gamma', type=float, default=0.95, help='The decay factor of learning rate')
parser.add_argument('--n_hidden', type=int, default=128, help='The Number of hidden neural used in MLP')
parser.add_argument('--K', type=int, default=1, help='Number of Gaussian kernel in GMM, default =1')
parser.add_argument('--K1', type=int, default=1, help='The Number of epoch for optimize MLP. Notes that we optimize MLP and W alternately. The default setting denotes to optimize MLP for one epoch then optimize W for two epochs.')
parser.add_argument('--K2', type=int, default=2, help='The Number of epoch for optimize W. Notes that we optimize MLP and W alternately. The default setting denotes to optimize MLP for one epoch then optimize W for two epochs.')
parser.add_argument('--save_name', type=str, default='/tmp')
opt = parser.parse_args()
if opt.task == 'non_celltype_GRN':
    if opt.setting == 'default':
        opt.beta = 1
        opt.alpha = 100
        opt.K1 = 1
        opt.K2 = 2
        opt.n_hidden = 128
        opt.gamma = 0.95
        opt.lr = 1e-4
        opt.lr_step_size = 0.99
        opt.batch_size = 64
    model = non_celltype_GRN_model(opt)
    if opt.setting == 'test':
        opt.beta = 1
        opt.alpha = 100
        opt.K1 = 1
        opt.K2 = 2
        opt.n_hidden = 128
        opt.gamma = 0.95
        opt.lr = 1e-4
        opt.lr_step_size = 0.99
        opt.batch_size = 64
        model = test_non_celltype_GRN_model(opt)
    model.train_model()
elif opt.task == 'celltype_GRN':
    if opt.setting == 'default':
        opt.beta = 0.01  # we found beta=0.01 alpha =1 perform better than alpha=0.1 beta=10 after paper submission
        opt.alpha = 1
        opt.K1 = 1
        opt.K2 = 2
        opt.n_hidden = 128
        opt.gamma = 0.95
        opt.lr = 1e-4
        opt.lr_step_size = 0.99
        opt.batch_size = 64
    model = celltype_GRN_model(opt)
    if opt.setting == 'test':
        opt.beta = 0.01
        opt.alpha = 1
        opt.K1 = 1
        opt.K2 = 2
        opt.n_hidden = 128
        opt.gamma = 0.95
        opt.lr = 1e-4
        opt.lr_step_size = 0.99
        opt.batch_size = 64
        model = test_celltype_GRN_model(opt)
    model.train_model()
elif opt.task == 'simulation':
    if opt.setting == 'default':
        opt.n_epochs = 120
        opt.beta = 1
        opt.alpha = 10
        opt.K1 = 1
        opt.K2 = 2
        opt.n_hidden = 128
        opt.gamma = 0.95
        opt.lr = 1e-4
        opt.lr_step_size = 0.99
        opt.batch_size = 64
    model = deepsem_generation(opt)
    model.train_model()
elif opt.task == 'embedding':
    if opt.setting == 'default':
        opt.n_epochs = 120
        opt.beta = 1
        opt.alpha = 10
        opt.K1 = 1
        opt.K2 = 2
        opt.n_hidden = 128
        opt.gamma = 0.95
        opt.lr = 1e-4
        opt.lr_step_size = 0.99
        opt.batch_size = 64
        opt.K = 1
    model = deepsem_embed(opt)
    model.train_model()
