# Popescu Ioana-Maria, 334CB

# (Current state, word)
Conf = (int, str)

# Get states from DFA that can reach a final state
def dfs(visited, graph, node):
    # for node in nodes:
    if node not in visited:
        visited.append(node)
        if node in graph:
            for neighbour in graph[node]:
                dfs(visited, graph, neighbour)
    # print("Result of dfs: ", visited)
    return visited
    # pass


# Reverse delta 
def reverse_delta(delta):
    rev = {}
    for key in delta.keys():
        for symbol, state in delta[key].items():
            if state not in rev:
                rev[state] = [key]
            else:
                rev[state].append(key)
    return rev


class DFA:
    # Initialize DFA fields
    def __init__(self, alphabet, token, initial_state, delta, sink, final_states):
        self.alphabet = alphabet
        self.token = token
        self.initialState: int = initial_state
        self.delta = delta
        self.sink = sink
        self.finalStates = final_states

    # Computes the next configuration from curr_config, if possible
    def step(self, conf: Conf) -> Conf:
        if conf[0] in self.delta and conf[1][0] in self.delta[conf[0]]:
            return self.delta[conf[0]][conf[1][0]], conf[1][1:]
        return -1, conf[1][1:]

    # Verifies if word is accepted by DFA
    def accept(self, word: str) -> bool:
        cnf = (self.initialState, word)
        while cnf[1]:
            if cnf[0] not in self.sink:
                return 0
            cnf = self.step(cnf)
        return cnf[0] in self.finalStates


# Splits the alphabet string in a list of each character
def alphabet_split(word):
    alphabet = []
    for i in range(0, len(word)):
        if word[i] == "\\" and word[i + 1] == "n":
            alphabet.append("\n")
            i += 1
        elif word[i - 1] == "\\" and word[i] == "n":
            continue
        else:
            alphabet.append(word[i])
    return alphabet


# Reads each field and returns the built dfa
def read_dfa(l):
    alphabet = alphabet_split(l[0][:-1])
    token = l[1][:-1]
    initial_state = int(l[2][:-1])

    delta = {}
    # Build Delta
    for i in range(3, len(l) - 1):
        transition = l[i].split(',')
        if int(transition[0]) in delta:
            d = delta[int(transition[0])]
            if transition[1][1] != "\\":
                d[transition[1][1]] = int(transition[2])
            else:
                d["\n"] = int(transition[2])
        else:
            if transition[1][1] != "\\":
                delta[int(transition[0])] = {transition[1][1]: int(transition[2])}
            else:
                delta[int(transition[0])] = {"\n": int(transition[2])}

    # Build final states
    last_line = l[-1].split()
    final_states = []
    for x in last_line:
        x = int(x)
        final_states.append(x)

    # Get all states that aren't sink states
    visited = []
    rev_delta = reverse_delta(delta)
    for i in final_states:
        not_sink = dfs(visited, rev_delta, i)

    return DFA(alphabet, token, initial_state, delta, not_sink, final_states)


# Searches for the longest prefix match
def longest_prefix(dfa_list, word):
    aux = word
    length = len(word)
    rez = []
    # Goes through the word and returns the the longest words accepted
    while aux:
        if dfas_accept(dfa_list, aux, rez):
            aux = word[length:]
            length = len(word)
        else:
            aux = aux[:-1]
            length -= 1
    return rez


# Goes through all the dfas and checks if/what dfa accepts it
def dfas_accept(dfa_list, word, rez):
    for dfa in dfa_list:
        if dfa.accept(word):
            if word != "\n":
                rez.append(dfa.token + " " + word)
            else:
                rez.append(dfa.token + " \\" + "n")
            return 1
    return 0


# Functia pentru checker
def runlexer(lexer, finput, foutput):
    i = 0
    elems = []
    dfas = []
    with open(lexer, "r") as lex:
        lines = lex.readlines()
        while i < len(lines) and lines[i][-1] == "\n":
            while i < len(lines) and lines[i] != "\n":
                elems.append(lines[i])
                i += 1
            dfa = read_dfa(elems)
            dfas.append(dfa)
            i += 1
            elems.clear()

    with open(finput, "r") as input:
        word = ''.join(input.readlines())

    rez = longest_prefix(dfas, word)
    with open(foutput, "w") as output:
        for line in rez:
            output.write(line + "\n")
    lex.close()
    input.close()
    output.close()
