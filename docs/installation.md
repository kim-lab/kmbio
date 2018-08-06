# Installation

## Stable release

To install `kmbio`, run this command in your terminal:

```bash
conda install kmbio -c kimlab
```

This is the preferred method to install `kmbio`, as it will always install the most recent stable release.

If you don't have [conda] installed, this [Python installation guide] can guide
you through the process.

[conda]: https://conda.io
[Python installation guide]: https://conda.io/docs/user-guide/install/index.html

## From sources

The sources for `kmbio` can be downloaded from the [GitLab repo].

You can either clone the public repository:

```bash
git clone git://gitlab.com/kimlab/kmbio
```

Or download the [tarball]:

```bash
curl -OL https://gitlab.com/kimlab/kmbio/repository/master/archive.tar
```

Once you have a copy of the source, you can install it with:

```bash
python setup.py install
```

[GitLab repo]: https://gitlab.com/kimlab/kmbio
[tarball]: https://gitlab.com/kimlab/kmbio/repository/master/archive.tar
