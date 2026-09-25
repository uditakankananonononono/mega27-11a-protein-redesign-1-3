"""One strictly protein-group held-out benchmark; no evaluation leakage."""
from pathlib import Path
import json
import numpy as np
import torch
from scipy.stats import spearmanr, pearsonr
from sklearn.model_selection import GroupShuffleSplit
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import Ridge
from .features import load_graphs, tabular_features
from .model import ResidueGNN, collate


def metrics(y, prediction):
    return dict(n=len(y), mae=float(np.mean(np.abs(y - prediction))), rmse=float(np.sqrt(np.mean((y-prediction)**2))),
                spearman=None if np.std(prediction) == 0 else float(spearmanr(y, prediction).statistic),
                pearson=None if np.std(prediction) == 0 else float(pearsonr(y, prediction).statistic))


def run(data: Path, structures: Path, output: Path, epochs: int = 80, seed: int = 27):
    torch.manual_seed(seed)
    np.random.seed(seed)
    torch.set_num_threads(2)
    graphs, audit = load_graphs(data, structures)
    if len(graphs) < 20:
        raise ValueError(f'Too few aligned structure examples: {audit}')
    y = np.array([g['y'] for g in graphs])
    groups = np.array([g['pdb'] for g in graphs])
    train, test = next(GroupShuffleSplit(n_splits=1, test_size=.25, random_state=seed).split(y, y, groups))
    assert not set(groups[train]) & set(groups[test])
    features = np.array([tabular_features(g) for g in graphs])
    ridge = make_pipeline(StandardScaler(), Ridge(alpha=100.)).fit(features[train], y[train])
    baseline = np.full(len(test), y[train].mean())
    ridge_pred = ridge.predict(features[test])
    model = ResidueGNN()
    opt = torch.optim.AdamW(model.parameters(), lr=.003, weight_decay=.01)
    xtrain, atrain, ytrain = collate([graphs[i] for i in train])
    xtest, atest, _ = collate([graphs[i] for i in test])
    for epoch in range(epochs):
        model.train()
        for batch in torch.randperm(len(train)).split(32):
            opt.zero_grad()
            loss = torch.nn.functional.smooth_l1_loss(model(xtrain[batch], atrain[batch]), ytrain[batch])
            loss.backward()
            opt.step()
    model.eval()
    with torch.no_grad():
        prediction = model(xtest, atest).numpy()
    raw = __import__('pandas').read_csv(data).set_index('name')
    published = {}
    for method in ['GeoDDG-Seq_dir', 'GeoDDG-3D_dir', 'PremPS_dir', 'DDGun_dir', 'ACDC-NN_dir', 'DDMut_dir']:
        observed = np.array([raw.loc[graphs[i]['name'], method] for i in test], dtype=float)
        keep = np.isfinite(observed)
        published[method] = metrics(y[test][keep], observed[keep])
    result = dict(seed=seed, epochs=epochs, audit=audit,
                  split=dict(train=len(train), test=len(test), train_pdbs=sorted(set(groups[train])), test_pdbs=sorted(set(groups[test]))),
                  metrics={name: metrics(y[test], pred) for name, pred in [('train_mean', baseline), ('ridge', ridge_pred), ('gnn', prediction)]} | published,
                  training=dict(model='ResidueGNN',architecture='42 features, 48 hidden units, 2 residual message-passing layers, masked mean plus site readout',optimizer='AdamW',learning_rate=0.003,weight_decay=0.01,loss='SmoothL1',batch_size=32,torch_version=torch.__version__),
                  predictions=[dict(name=graphs[i]['name'], pdb=groups[i], experimental_ddg=float(y[i]), train_mean=float(baseline[j]), ridge=float(ridge_pred[j]), gnn=float(prediction[j])) for j,i in enumerate(test)])
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(result, indent=2, allow_nan=False) + '\n')
    return result
