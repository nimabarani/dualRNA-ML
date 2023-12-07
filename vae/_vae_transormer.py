import torch
from torch.optim import Adam
from sklearn.base import BaseEstimator, TransformerMixin
from torch.utils.data import DataLoader

from ._vae import VAE
from ._data_set import CustomDataset
from ._vae_loss import VAELoss


class VAETransformer(BaseEstimator, TransformerMixin):
    def __init__(
        self,
        input_dim=20,
        latent_dim=50,
        num_epochs=200,
        batch_size=256,
        alpha=1,
        learning_rate=0.002,
        verbose=None,
    ):
        self.__device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
        self.input_dim = input_dim
        self.latent_dim = latent_dim
        self.num_epochs = num_epochs
        self.batch_size = batch_size
        self.alpha = alpha
        self.learning_rate = learning_rate
        self.verbose = verbose

    def fit(self, X, y=None):
        # Convert X to a torch tensor and move it to the appropriate device
        self._vae = VAE(self.input_dim, self.latent_dim).to(self.__device)
        X = torch.tensor(X, dtype=torch.float32).to(self.__device)
        dataset = CustomDataset(X)
        data_loader = DataLoader(dataset, batch_size=self.batch_size, shuffle=True)

        optimizer = Adam(self._vae.parameters(), lr=self.learning_rate)
        vae_loss = VAELoss()
        for epoch in range(self.num_epochs):
            for batch in data_loader:
                optimizer.zero_grad()
                recon_batch, z_mean, z_log_var = self._vae(batch)
                loss = vae_loss(
                    batch,
                    recon_batch,
                    z_mean,
                    z_log_var,
                    self.input_dim,
                    self.alpha,
                )
                loss.backward()
                optimizer.step()

            if self.verbose:
                print(f"Epoch [{epoch+1}/{self.num_epochs}], Loss: {loss.item()}")
        return self

    def transform(self, X, y=None):
        # Convert X to a torch tensor and move it to the appropriate device
        X = torch.tensor(X, dtype=torch.float32).to(self.__device)

        self._vae.eval()
        with torch.no_grad():
            latent_params = self._vae.encoder(X).detach()
            mu, _ = latent_params.chunk(2, dim=1)
        return mu.cpu().numpy()
