def special(num_genes, total_gene_list):
     fin_genes = []
     while len(fin_genes) != num_genes:
         gene_idx = np.random.randint(0, len(total_gene_list), 1)[0]
         if b'.' not in total_gene_list[gene_idx] and total_gene_list[gene_idx] not in bad_genes and total_gene_list[
 gene_idx] not in fin_genes:
             a = total_gene_list[gene_idx]
             #print(a)
             a =a.decode("utf-8")
             fin_genes.append(a)
     return fin_genes

      def clean(a):
          a = a.replace(" ","")
          a = a.replace("\t","")
          a = a.split(",")
          return a
