import torch
import torch.optim as optim
from sklearn.base import BaseEstimator, TransformerMixin
from torch import nn
from torch.utils.data import DataLoader

from ._vae import VAE
from ._data_set import CustomDataset


class VAETransformer(BaseEstimator, TransformerMixin):
    def __init__(
        self,
        input_dim,
        latent_dim=50,
        num_epochs=100,
        batch_size=128,
        verbose=None,
    ):
        self.__device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

        self.vae = VAE(input_dim, latent_dim).to(self.__device)
        self.input_dim = input_dim
        self.latent_dim = latent_dim
        self.num_epochs = num_epochs
        self.batch_size = batch_size
        self.verbose = verbose
        

    def fit(self, X, y=None):
        # Convert X to a torch tensor and move it to the appropriate device
        # X = torch.tensor(X, dtype=torch.float32).to(self.device)
        # input_dim = X.shape[1]

        dataset = CustomDataset(X)
        data_loader = DataLoader(dataset, batch_size=self.batch_size, shuffle=True)

        self.vae = VAE(self.input_dim, self.latent_dim).to(self.__device)

        optimizer = optim.Adam(self.vae.parameters(), lr=0.001)

        # Assuming data_loader is already defined
        for epoch in range(self.num_epochs):
            for i, batch in enumerate(data_loader):
                # Prepare inputs
                # Assuming each batch consists of input data and corresponding labels
                inputs = batch
                # Flatten the input if needed
                inputs = inputs.view(inputs.size(0), -1).to(self.__device)

                # Forward pass
                reconstructed_inputs, mu, logvar = self.vae(inputs)

                # Compute the reconstruction loss (binary cross-entropy)
                reconstruction_loss = nn.functional.binary_cross_entropy(
                    reconstructed_inputs, inputs, reduction="sum"
                )

                # Compute the KL divergence loss
                kl_divergence_loss = -0.5 * torch.sum(
                    1 + logvar - mu.pow(2) - logvar.exp()
                )

                # Compute the total loss
                total_loss = reconstruction_loss + kl_divergence_loss

                # Backward pass and optimization
                optimizer.zero_grad()
                total_loss.backward()
                optimizer.step()

            # Print loss for each epoch
            if self.verbose:
                print(f"Epoch [{epoch+1}/{self.num_epochs}], Loss: {total_loss.item()}")

        return self

    def transform(self, X, y=None):
        # Convert X to a torch tensor and move it to the appropriate device
        X = torch.tensor(X, dtype=torch.float32).to(self.device)
        self.vae.eval()

        with torch.no_grad():
            # Pass through the encoder and detach to remove it from computation graph
            latent_params = self.vae.encoder(X).detach()
            mu, _ = latent_params.chunk(2, dim=1)
        return mu.cpu().numpy()
