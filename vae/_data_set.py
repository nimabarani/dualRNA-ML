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
