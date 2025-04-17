#TODO import

#TODO need to specify localcolabfold directory in the container bind mounts
class Simulation:
    def __init__(self, protein):
        self.protein = protein

    def cleanup(self):
        pass

    #TODO openmm with a warning
    def simulate(self, use_bioemu=True, relax=False):
        pass

    def save(self):
        pass

    
