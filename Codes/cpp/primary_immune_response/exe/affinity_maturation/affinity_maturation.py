import random

class SequenceAgent:
    def __init__(self, unique_id, sequence):
        self.unique_id = unique_id
        self.sequence = sequence

    def step(self, agents, birth_rate, death_rate):
        # Calculate the total birth and death rates
        total_birth_rate = birth_rate * len(agents)
        total_death_rate = death_rate * len(agents)

        # Calculate the total rate of all possible events
        total_rate = total_birth_rate + total_death_rate

        # Determine the time to the next event
        time_to_next_event = random.expovariate(total_rate)

        # Determine the type of event (birth or death) based on the probabilities
        if random.uniform(0, total_rate) < total_birth_rate:
            self.birth(agents)
        else:
            self.death(agents)

        return time_to_next_event

    def birth(self, agents):
        # Randomly select a neighbor
        neighbor = random.choice(agents)

        # Randomly select a position within the sequence
        position = random.randint(0, len(self.sequence) - 1)

        # Replace the element at the selected position with the neighbor's value
        mutated_sequence = list(self.sequence)
        mutated_sequence[position] = neighbor.sequence[position]

        # Create a new agent with the mutated sequence
        child = SequenceAgent(len(agents), mutated_sequence)
        agents.append(child)

    def death(self, agents):
        agents.remove(self)

# Parameters
N = 10  # Number of agents
sequence_length = 5  # Length of each sequence
num_steps = 10  # Number of simulation steps
birth_rate = 0.1  # Rate of birth events
death_rate = 0.2  # Rate of death events

# Create initial agents with random sequences
agents = []
for i in range(N):
    sequence = [random.randint(0, 19) for _ in range(sequence_length)]
    agent = SequenceAgent(i, sequence)
    agents.append(agent)

# Simulation loop using Gillespie algorithm
time = 0
while time < num_steps:
    # Calculate the total birth and death rates
    total_birth_rate = birth_rate * len(agents)
    total_death_rate = death_rate * len(agents)

    # Calculate the total rate of all possible events
    total_rate = total_birth_rate + total_death_rate

    # Determine the time to the next event
    time_to_next_event = random.expovariate(total_rate)

    # Check if the next event occurs within the simulation time
    if time + time_to_next_event > num_steps:
        break

    # Increment the time
    time += time_to_next_event

    # Determine the type of event (birth or death) based on the probabilities
    if random.uniform(0, total_rate) < total_birth_rate:
        chosen_agent = random.choice(agents)
        chosen_agent.birth(agents)
    else:
        chosen_agent = random.choice(agents)
        chosen_agent.death(agents)

    print("Time:", time)
    for agent in agents:
        print("Agent", agent.unique_id, "Sequence:", agent.sequence)
    print()
