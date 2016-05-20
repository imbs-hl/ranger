## ranger python bindings

### Overview

These are pretty minimal python bindings for ranger - they work but don't expose much functionality. Essentially only a RandomForestClassifier class with sklearn-like API is exposed.


### Installation

```
./waf configure
./waf
./waf install
```

### Example usage

```python
import pyranger
import numpy as np
from sklearn import datasets

iris = datasets.load_iris()
X = iris.data
y = iris.target

rf = pyranger.RandomForestClassifier()
rf.fit(X, y)
y_pred = rf.predict(X)

print(np.where(y_pred != y))
```