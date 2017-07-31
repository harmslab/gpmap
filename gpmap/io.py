import json
import pickle
import pandas as pd

class GPMIO(object):

    @classmethod
    def read_excel(cls, fname, **kwargs):
        """"""
        pd.read_excep
    @classmethod
    def read_csv(cls, fname, **kwargs):
        """"""


    # ----------------------------------------------------------
    # Class method to load from source
    # ----------------------------------------------------------

    @classmethod
    def read_json(cls, filename, **kwargs):
        """Load a genotype-phenotype map directly from a json file.
        The JSON metadata must include the following attributes

        Note
        ----
        Keyword arguments override input that is loaded from the JSON file.
        """
        # Open, json load, and close a json file
        f = open(filename, "r")
        data = json.load(f)
        f.close()

        # Grab all properties from data-structure
        necessary_args = ["wildtype", "genotypes", "phenotypes"]
        options = {
            "genotypes" : [],
            "phenotypes" : [],
            "wildtype" : [],
            "stdeviations": None,
            "log_transform": False,
            "mutations": None,
            "n_replicates": 1,
            "logbase": np.log10
        }
        # Get all options for map and order them
        for key in options:
            # See if options are in json data
            try:
                options[key] = data[key]
            except:
                pass
        # Override any properties with manually entered kwargs passed directly into method
        options.update(kwargs)
        args = []
        for arg in necessary_args:
            val = options.pop(arg)
            args.append(val)
        # Create an instance
        gpm = cls(args[0], args[1], args[2], **options)
        return gpm

    def to_json(self, filename):
        """Write genotype-phenotype map to json file.
        """
        # Get metadata and make sure values are lists.
        data = self.metadata
        for key, val in data.items():
            if hasattr(val, "__iter__") and type(val) != dict:
                data[key] = list(val)
        # Write to file
        with open(filename, "w") as f:
            json.dump(data, f)

    def to_csv(self, fname, items, sep="\t", **kwargs):
        """Write out items to a tab-separated file.
        """
        # write a mapping dictionary to file
        if type(items) is dict:
            with open(fname, "w") as f:
                for key, value in items.items():
                    f.write(str(key) + str(sep) + str(value) + "\n")
        # write a list of items to file.
        else:
            with open(fname, "w") as f:
                nrows = len(items[0])
                ncols = len(items)
                for i in range(nrows):
                    row = []
                    for j in range(ncols):
                        row.append(items[j][i])
                    f.write(sep.join(row))
