import sys


def merge_two_dicts(x, y):
    z = x.copy()
    z.update(y)
    return z


def str_to_list(key):
    key_to_list = key[1:-1].split(", ")
    for i in range(0, len(key_to_list)):
        key_to_list[i] = int(key_to_list[i])
    return key_to_list


def find_index(l, check, sink):
    for el in check:
        if str(l) == el[0]:
            return el[2]
    return sink


class DFA:
    # Initialize DFA fields
    def __init__(self, alphabet, initial_state, delta, nr_states, final_states):
        self.alphabet = alphabet
        self.initialState: int = initial_state
        self.delta = delta
        self.nr_states = nr_states
        self.finalStates = final_states


class NFA:
    # Initialize NFA fields
    def __init__(self, alphabet, nr_states, initial_state, delta, final_state):
        self.alphabet = alphabet
        self.nr_states = nr_states
        self.initial_state = initial_state
        self.delta = delta
        self.final_state = final_state


class Expr:
    name: str

    def __init__(self, name):
        self.name = name

    def nr_params(self):
        if self.name == "UNION" or self.name == "CONCAT":
            return 2
        elif self.name == "STAR" or self.name == "PLUS":
            return 1
        else:
            return 0


class Symbol(Expr):
    char: str

    def __init__(self, c):
        self.char = c

    def symbolNFA(self, c, all_states):
        alphabet = [c]
        nr_states = 2
        initial_state = all_states
        delta = {initial_state: {c: [initial_state + 1]}}
        final_state = initial_state + 1
        return NFA(alphabet, nr_states, initial_state, delta, final_state)


class Star(Expr):

    def __init__(self, expr1):
        self.expr1 = expr1

    def starNFA(self, nfa1: NFA, all_states: int):
        alphabet = nfa1.alphabet
        nr_states = nfa1.nr_states + 2
        initial_state = all_states
        final_state = all_states + 1

        delta = nfa1.delta

        delta[initial_state] = {"": [nfa1.initial_state, final_state]}
        delta[nfa1.final_state] = {"": [nfa1.initial_state, final_state]}

        return NFA(alphabet, nr_states, initial_state, delta, final_state)


class Plus(Expr):
    def __init__(self, expr1):
        self.expr1 = expr1

    def plusNFA(self, nfa1: NFA, all_states: int):
        alphabet = nfa1.alphabet
        nr_states = nfa1.nr_states + 2
        initial_state = all_states
        final_state = all_states + 1

        delta = nfa1.delta

        delta[initial_state] = {"": [nfa1.initial_state]}
        delta[nfa1.final_state] = {"": [nfa1.initial_state, final_state]}

        return NFA(alphabet, nr_states, initial_state, delta, final_state)


class Concat(Expr):

    def __init__(self, expr1, expr2):
        self.expr1 = expr1
        self.expr2 = expr2

    def concNFA(self, nfa1: NFA, nfa2: NFA):
        alphabet = nfa1.alphabet
        alphabet.extend(nfa2.alphabet)
        nr_states = nfa1.nr_states + nfa2.nr_states
        initial_state = nfa1.initial_state
        delta = merge_two_dicts(nfa1.delta, nfa2.delta)
        delta[nfa1.final_state] = {"": [nfa2.initial_state]}
        final_state = nfa2.final_state
        return NFA(alphabet, nr_states, initial_state, delta, final_state)


class Union(Expr):

    def __init__(self, expr1, expr2):
        self.expr1 = expr1
        self.expr2 = expr2

    def unionNFA(self, nfa1: NFA, nfa2: NFA, all_states: int):
        alphabet = nfa1.alphabet
        alphabet.extend(nfa2.alphabet)
        nr_states = nfa1.nr_states + nfa2.nr_states + 2
        initial_state = all_states
        final_state = all_states + 1

        delta = merge_two_dicts(nfa1.delta, nfa2.delta)
        delta[initial_state] = {"": [nfa1.initial_state, nfa2.initial_state]}

        delta[nfa1.final_state] = {"": [final_state]}
        delta[nfa2.final_state] = {"": [final_state]}

        return NFA(alphabet, nr_states, initial_state, delta, final_state)


# DFS for computing epsilon-closures
def dfs(visited, graph, node):
    if node not in visited:
        visited.append(node)
        if node in graph:
            for item in graph[node].items():
                if item[0] == '':
                    for elem in item[1]:
                        dfs(visited, graph, elem)


# Get next state from epsilon closure
def next_states(nfa, e_closures, ll):
    result = []
    final = []
    alp = nfa.alphabet
    for i in range(0, len(nfa.alphabet)):
        l = []
        next = []
        for elem in ll:
            for nr in e_closures[elem]:
                if nr in nfa.delta and alp[i] in nfa.delta[nr]:
                    aux = nfa.delta[nr]
                    aux1 = aux[alp[i]]
                    l.extend(e_closures[aux1[0]])
                    next.extend(aux1)
        next = list(dict.fromkeys(next))
        next.sort()
        final.append(next)
        l = list(dict.fromkeys(l))
        l.sort()
        result.append(l)
    return (final, result)


# Get all dfa transitions
def dfs_transitions(visited, delta, group, e_closures, nfa, result):
    if group not in visited:
        visited.append(group)
        dest = next_states(nfa, e_closures, group)
        state = []
        for elems in group:
            state.extend(e_closures[elems])
        state = list(dict.fromkeys(state))
        state.sort()
        result[str(state)] = dict(zip(nfa.alphabet, dest[1]))
        for next_list in dest[0]:
            dfs_transitions(visited, delta, next_list, e_closures, nfa, result)


def regex_dfa(regex_file, dfa_file):

    with open(regex_file, "r") as rgx:
        inp = rgx.readlines()[0]
    rgx.close()
    all_states = 0
    stack = []
    nfas = []
    line = inp.split()

    # Regex Stack
    for token in line:
        stack.append(token)

    while stack:
        elem = stack.pop()
        expr = Expr(elem)

        # Symbol
        if expr.nr_params() == 0:
            symbol = Symbol(elem)
            nfas.append(symbol.symbolNFA(elem, all_states))
            all_states += 2
        # Star & Plus
        elif expr.nr_params() == 1:
            if elem == "STAR":
                star = Star(nfas.pop())
                nfas.append(star.starNFA(star.expr1, all_states))
                all_states += 2
            else:
                plus = Plus(nfas.pop())
                nfas.append(plus.plusNFA(plus.expr1, all_states))
                all_states += 2
        # Concat & Union
        else:
            if elem == "CONCAT":
                conc = Concat(nfas.pop(), nfas.pop())
                nfas.append(conc.concNFA(conc.expr1, conc.expr2))
            else:
                union = Union(nfas.pop(), nfas.pop())
                nfas.append(union.unionNFA(union.expr1, union.expr2, all_states))
                all_states += 2

    # Remove duplicates from alphabet
    for i in range(0, len(nfas)):
        nfas[i].alphabet = list(dict.fromkeys(nfas[i].alphabet))

    # Compute Epsilon Closures
    epsilon_closures = []
    for i in range(0, all_states):
        visited = []
        dfs(visited, nfas[0].delta, i)
        visited.sort()
        epsilon_closures.append(visited)
    alp = nfas[0].alphabet

    nfa = nfas[0]
    visited = []
    result = {}
    dfs_transitions(visited, nfa.delta, [nfas[0].initial_state], epsilon_closures, nfa, result)

    index = 0
    check = []
    for key in result:
        if key != '[]':
            l = str_to_list(key)
            check.append((key, l, index))
        else:
            check.append((key, [], index))
        index += 1


    final_delta = {}
    for key in result.items():
        trans = {}
        for x in key[1].items():
            trans[x[0]] = find_index(x[1], check, index - 1)

        final_delta[find_index(key[0], check, index - 1)] = trans

    # Compute final states
    dfa_finalstates = []
    for elem in check:
        if nfas[0].final_state in elem[1]:
            dfa_finalstates.append(elem[2])

    dfa = DFA(alp, 0, final_delta, index, dfa_finalstates)

    # Write output
    with open(dfa_file, "w") as output:
        for letter in dfa.alphabet:
            output.write(letter)
        output.write("\n" + str(dfa.nr_states) + "\n" + str(dfa.initialState) + "\n")

        for i in range(0, len(dfa.finalStates)):
            if i == len(dfa.finalStates) - 1:
                output.write(str(dfa.finalStates[i]))
            else:
                output.write(str(dfa.finalStates[i]) + " ")
        output.write("\n")
        for x in dfa.delta.items():
            for y in x[1].items():
                output.write(str(x[0]) + ",'" + y[0] + "'," + str(y[1]) + "\n")
    output.close()


if __name__ == '__main__':
    regex_dfa(sys.argv[1], sys.argv[2])
