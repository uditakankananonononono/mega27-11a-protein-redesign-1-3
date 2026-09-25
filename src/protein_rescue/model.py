"""Small residue graph network for experimentally grounded stability prediction."""
import torch
from torch import nn


class ResidueGNN(nn.Module):
    def __init__(self, input_dim=42, hidden=48):
        super().__init__()
        self.embedding = nn.Linear(input_dim, hidden)
        self.conv1 = nn.Linear(hidden, hidden)
        self.conv2 = nn.Linear(hidden, hidden)
        self.out = nn.Sequential(nn.Linear(hidden * 2, hidden), nn.ReLU(), nn.Linear(hidden, 1))

    def forward(self, x, adjacency):
        h = torch.relu(self.embedding(x))
        h = torch.relu(self.conv1(adjacency @ h)) + h
        h = torch.relu(self.conv2(adjacency @ h)) + h
        pooled = torch.cat([h.mean(dim=1), (h * x[:, :, 40:41]).sum(dim=1)], dim=1)
        return self.out(pooled).flatten()


def collate(graphs):
    n = max(len(g['x']) for g in graphs)
    x = torch.zeros((len(graphs), n, 42), dtype=torch.float32)
    a = torch.zeros((len(graphs), n, n), dtype=torch.float32)
    for i, g in enumerate(graphs):
        k = len(g['x'])
        x[i, :k] = torch.from_numpy(g['x'])
        a[i, :k, :k] = torch.from_numpy(g['a'])
    y = torch.tensor([g['y'] for g in graphs], dtype=torch.float32)
    return x, a, y
