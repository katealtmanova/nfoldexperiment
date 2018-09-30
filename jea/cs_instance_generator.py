import random, copy
from collections import defaultdict

def gen_strings(n,k,r,Sigma):
    # n = length of strings
    # k = number of strings
    # Sigma  = alphabet (list)
    # r = ratio n/\alpha, so we want to make \alpha changes
    alpha = int(n/r)
    
    master_string = [random.choice(Sigma) for i in range(n)]
    strings = [copy.copy(master_string) for i in range(k)]
    
    for i in range(alpha):
        kp = random.choice(range(k))
        pos = random.choice(range(n))
        Sig = copy.copy(Sigma)
        old_char = strings[kp][pos]
        Sig.remove(old_char)
        char = random.choice(Sig)
        strings[kp][pos] = char
    
    strings = ["".join(s) for s in strings]
    
    #out_lines = []
    #out_lines.append(str(len(Sigma)))
    #out_lines.append(str(k))
    #out_lines.append(str(n))
    #out_lines.extend(Sigma)
    #out_lines.extend(strings)
    
    out = {}
    out["Sigma_len"] = len(Sigma)
    out["k"] = k
    out["n"] = n
    out["Sigma"] = Sigma
    out["strings"] = strings
    
    return out

    
    

class Preprocess_data_to_dictionary:  
    
    def __init__(self,data):
        
        self.data = data
        #read the header
        self.alphabet_size = int(data["Sigma_len"])
        self.num_of_strings = int(data["k"])
        self.len_of_strings = int(data["n"])
        self.make_array_of_data()
        self.create_dictionary_from_array()
    
    def return_dictionary(self):
        
        #make an array of the strings from the file
        self.make_array_of_data()
        
        #create the dictionary from the saved array of the strings
        self.create_dictionary_from_array()
     
        return self.dict
        
    def make_array_of_data(self):
        #make list from the characters of an alphabet
        self.alphabet_list = self.data["Sigma"]
        
        self.data_array = self.data["strings"]
        #read the data --> to an array
        #self.data_array = [[0 for x in range(self.len_of_strings)] for y in range(self.num_of_strings)]         
        
        #for i in range(len(self.data["strings"])):
        #    for j in range(self.len_of_strings):
        #        self.data_array[i][j] = self.data["strings"][i][j]
         
    def create_dictionary_from_array(self):
        self.d = {}
        
        #for each columns
        for i in range(self.len_of_strings):
            
            s = ""        
            
            #build a string of a current column
            for j in range(self.num_of_strings):
                s += self.data_array[j][i]
                            
            #check if the string is already a key in the dictionary
            if s in self.d:
                self.d[s] += 1 #increment val
            else:
                self.d[s] = 1 #add key with value = 1
        
        # Preprocessing - reduce size of alphabet
        #print "old |\Sigma|:", self.alphabet_size
        #print "old #cols:", len(self.dict)
        d2 = defaultdict(int)
        for col in self.d:
            chars = list(set(col))
            new_chars = list(range(len(chars)))
            mapping = {}
            for (i,c) in enumerate(chars):
                mapping[c] = new_chars[i]
            new_col = "".join((str(mapping[c]) for c in col))
            d2[new_col] += self.d[col]
        complete_new_chars = set()
        for new_col in d2:
            complete_new_chars.update(new_col)
        complete_new_chars = list(complete_new_chars)
        
        self.d = d2
        self.alphabet_list = complete_new_chars
        self.alphabet_size = len(complete_new_chars)
        #print "new |\Sigma|:", self.alphabet_size
        #print "new #cols:", len(self.dict)


                
def getNFoldIPinstance(p, distance):
    from sage.matrix.constructor import matrix
    from sage.modules.free_module_element import vector
    # Basic preprocessing - delete obvious columns
    for e in p.alphabet_list:
        if e*p.num_of_strings in p.d:
            p.d.pop(e*p.num_of_strings) 
    
    B = sum(v for v in p.d.values())
    str_cols = p.d.keys()
    cols = len(str_cols)
    D_cols = []
    u = []
    # normal columns
    for col in str_cols:
        for e in p.alphabet_list:
            mask = tuple((1-int(c == e) for c in col))
            D_cols.append(mask)
            
    # for blanks
    D_cols.append((0,)*p.num_of_strings)
    t_noslacks = len(D_cols)
    # for slacks in D
    for i in range(p.num_of_strings):
        D_cols.append((0,)*i + (1,) + (0,)*(p.num_of_strings - i -1))
    t = len(D_cols)
    D_t = matrix(D_cols)
    D = D_t.transpose()
    A = matrix((1,)*t_noslacks + (0,)*p.num_of_strings)
    n = cols
    #distance = int(n/r + distance_factor*(k-1)*(n/r))
    distance = min(B, distance)
    l = [vector((0,)*t)]*n
    w = [vector((0,)*(p.alphabet_size*cols) + (1,) + (0,)*p.num_of_strings)]*n
    u = []
    for i in range(n):
        if i < n-1:
            bru = (0,)*(i*p.alphabet_size) + (B,)*p.alphabet_size + (0,)*((cols-i-1)*p.alphabet_size) + (B,) + (0,)*p.num_of_strings
        else:
            bru = (0,)*(i*p.alphabet_size) + (B,)*p.alphabet_size + (0,)*((cols-i-1)*p.alphabet_size) + (B,) + (B,)*p.num_of_strings
        u.append(vector(bru))
    x = [vector((0,)*(p.alphabet_size*cols) + (int(p.d[col]),) + (int(i == len(str_cols)-1)*distance,)*p.num_of_strings) for (i,col) in enumerate(str_cols)]
    b = [vector((distance,)*p.num_of_strings)] + [vector((p.d[col],)) for col in str_cols]
    
    
    return A, D, n, b, l, u, w, x

def gen_cs_instance(n,k,r,sigma,distance_factor):
    # real distance is between 0 and int(n/r)
    # because we made int(n/r) changes (across all strings) in total
    # distance_factor is between 0 and 1
    # so distance = int(distance_factor*(n/r))
    Sigma = [str(char) for char in range(sigma)]
    out = gen_strings(n,k,r,Sigma)
    dict_data = Preprocess_data_to_dictionary(out)
    distance = int(distance_factor*(n/r))
    inst = getNFoldIPinstance(dict_data, distance)
    return inst
    
