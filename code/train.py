import argparse
import os
import os.path as osp
import pytorch_lightning as pl
import lightning_lite
import torch.nn.functional as F
from utils import utils
from data.data import get_dl, gety
from models.GTCN import GTCN
from pytorch_lightning.loggers import WandbLogger
from pytorch_lightning.callbacks import TQDMProgressBar
import torch
import wandb


os.environ['CUDA_VISIBLE_DEVICES'] = '0'

parser = argparse.ArgumentParser(description='Traing GTCN')
parser.add_argument('--gpu', '-g', dest='gpu', default=0)
parser.add_argument('--seed', '-s', dest='seed', default=124, type=int)
parser.add_argument('--fold', '-f', dest='fold', default=0, type=int)
args = parser.parse_args()

lightning_lite.utilities.seed.seed_everything(seed=args.seed, workers=True)

yaml_path = osp.join(os.getcwd(), 'test.yaml')
config = utils.load_yaml(yaml_path)

model = GTCN(config)



def train(config, data_path, Total_num, splitpara, logger=True):

    print(config)

    train_dl, valid_dl = get_dl(config, data_path, Total_num, splitpara)
    print(f'train: {len(train_dl)} valid: {len(valid_dl)}')

    print('Start training:')
    EPOCHS = config.epochs
    checkpoint_callback = pl.callbacks.ModelCheckpoint(monitor='valid_mae', mode='min')
    pcb = TQDMProgressBar(refresh_rate=1)
    trainer = pl.Trainer(accelerator='gpu', devices=1,
                     strategy=None,
                     max_epochs=EPOCHS, 
                     callbacks=[checkpoint_callback, pcb],
                     logger=logger,
                     precision=16,
                     gradient_clip_val=config.gradient_clip_val,
                     gradient_clip_algorithm="value"
                    )

    trainer.fit(model, train_dataloaders=train_dl, val_dataloaders=valid_dl)

    print('Predict valid:')
    y_p = trainer.predict(model, valid_dl)
    y_p = torch.cat(y_p, dim=0)
    y = gety(valid_dl)

    score1 = F.l1_loss(y_p, y).item()

    msgs = {'VALID MAE':score1}
    print(f'Valid MAE: {score1:.4f}')
    print("weights saved:", trainer.log_dir)
    return msgs


if __name__ == "__main__":
    quick_run = False                                   
    
    dataset_path = r'F:\data'
    train(config=config, data_path=dataset_path, Total_num=900, splitpara=[720, 180])

