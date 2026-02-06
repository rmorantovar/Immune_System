import sys
from functions_1 import*

# ------- CLASSES -------
class Clone:
    def __init__(self, id, sequence, e, f, parent, counts=None):
        self.id = id
        self.sequence = sequence
        self.counts = counts if counts is not None else [1]
        self.e = e
        self.f = f
        self.parent = parent

class Population:
    def __init__(self, N_clones, N_cells, l, b, d, m, motif, pop0 = False, sequences_pop0 = [], ns_pop0 = []):
        self.N_clones = [N_clones]
        self.N_cells = [N_cells]
        self.l = l
        self.b = b
        self.d = d
        self.m = m
        self.clones = []
        self.next_id = N_clones+1
        self.time = 0
        self.motif = motif
        self.pop0 = pop0
        self.sequences_pop0 = sequences_pop0
        self.ns_pop0 = ns_pop0
        self.create_initial_population()

    def create_initial_population(self):
        if self.pop0:
            for n in range(len(self.sequences_pop0)):
                seq = self.sequences_pop0[n]
                counts = [self.ns_pop0[n]]
                e =  self.energy(seq)
                self.clones.append(Clone(id=n+1, sequence=seq, e = e, f = self.fitness(e = e), parent = 0, counts = counts))
        else:    
            for n in range(self.N_clones[0]):
                seq = np.random.randint(0, 20, self.l)
                e =  self.energy(seq)
                self.clones.append(Clone(id=n+1, sequence=seq, e = e, f = self.fitness(e = e), parent = 0, counts = counts))
            
    def fitness(self, e):
        return  1/(1+np.exp(e+14))

    def energy(self, seq):
        return calculate_energy(self.motif, seq) - 24

    def update_energies(self):
        for n in range(self.next_id-1):
            self.clones[n].e = calculate_energy(self.motif, self.clones[n].sequence) - 24

    def step(self):
        rates = []
        for clone in self.clones:
            # Birth, death, and mutation rates for each individual of the clone
            rates.append(self.b * clone.counts[-1]*(1 - self.N_cells[-1]/1e4) * clone.f)
            rates.append(self.d * clone.counts[-1] * (1 - clone.f))
            rates.append(self.m * clone.counts[-1])

        total_rate = sum(rates)
        dt = -np.log(np.random.rand()) / total_rate
        self.time += dt

        event = np.searchsorted(np.cumsum(rates)/total_rate, np.random.rand())
        #print(event)
        selected_clone_index = event // 3
        other_clones_index = [i for i in range(self.next_id-1)]
        other_clones_index.remove(selected_clone_index)
        selected_event = event % 3
        selected_clone = self.clones[selected_clone_index]
        # Birth event
        if selected_event == 0:
            # print('birth!')
            selected_clone.counts.append(selected_clone.counts[-1] +  1)
            self.N_cells.append(self.N_cells[-1] + 1)
        # Death event
        elif selected_event == 1:
            # print('death!')
            selected_clone.counts.append(selected_clone.counts[-1] -  1)
            self.N_cells.append(self.N_cells[-1] - 1)
        # Mutation event
        else:
            # print('mutation!')
            mutated_sequence = selected_clone.sequence.copy()
            mutated_sequence[np.random.randint(0, self.l)] = np.random.randint(0, 20)
            prev_counts = [0 for _ in selected_clone.counts]
            prev_counts.append(1)
            e_mut = self.energy(mutated_sequence)
            self.clones.append(Clone(id=self.next_id, sequence=mutated_sequence, e = e_mut, f = self.fitness(e_mut), parent = selected_clone.id, counts = prev_counts))
            self.next_id += 1
            self.N_clones.append(self.N_clones[-1] + 1)
            self.N_cells.append(self.N_cells[-1] + 1)
            selected_clone.counts.append(selected_clone.counts[-1])

        for j in other_clones_index:
            self.clones[j].counts.append(self.clones[j].counts[-1])
        return total_rate


    def simulate(self, T):        
        states = []
        times = [self.time]
        total_rate = 1
        while ( (self.time < T) and  (self.N_cells[-1]!=0)):
            total_rate = self.step()
            times.append(self.time)

        return self.clones, self.N_clones, self.N_cells, times, self.next_id


