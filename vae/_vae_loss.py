import torch
from torch import nn
import torch.nn.functional as F


class VAELoss(nn.Module):
    def __init__(self):
        super().__init__()

    def forward(self, x, x_decoded, z_mean, z_log_var, original_dim, alpha):
        reconstruction_loss = original_dim * F.l1_loss(x, x_decoded, reduction="sum")
        kl_loss = -0.5 * torch.sum(1 + z_log_var - z_mean.pow(2) - z_log_var.exp())
        return reconstruction_loss + alpha * kl_loss
