import torch
from torch.utils.data import Dataset


class CustomDataset(Dataset):
    def __init__(self, data, transform=None):
        self.data = data
        self.transform = transform

    def __len__(self):
        return len(self.data)

    def __getitem__(self, index):
        features = torch.tensor(self.data[index], dtype=torch.float32)

        if self.transform:
            features = self.transform(features)

        return features


# class VAELoss(nn.Module):
#     def forward(self, x, x_recon, mu, logvar):
#         recon_loss = nn.functional.binary_cross_entropy(x_recon, x, reduction="sum")
#         kld_loss = -0.5 * torch.sum(1 + logvar - mu.pow(2) - logvar.exp())
#         return recon_loss + kld_loss
