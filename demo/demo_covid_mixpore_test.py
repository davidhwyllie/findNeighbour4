""" illustrates use of findNeighbour4 with covid samples, doing mixpore test
assumes a findNeighbour4 server is running
"""

if __name__ == '__main__':
    
    import random
    from fn4client import fn4Client

    # instantiate client
    fn4c = fn4Client("http://findneighbours04.unix.phe.gov.uk:5023")      # expects operation on local host; pass baseurl if somewhere else.

    existing_guids = set(fn4c.guids())
    print("There are {0} existing guids".format(len(existing_guids)))

   
    # get the neighbours of my_special_samples
    for this_guid in existing_guids:

        res = fn4c.guid2neighbours(this_guid, threshold = 1)        # find neighbours within 1 SNV
        related_samples = set()
        for related_sample, distance in res:
            related_samples.add(related_sample)
        print(this_guid, "has {0} neighbours".format(len(related_samples)))

        # if there are more than 20 related samples, randomly sample
        if len(related_samples)  > 20:
            related_samples = set(random.sample(related_samples, 20))
        related_samples.add(this_guid)
            
        # note that the for_msa call can return
       # print("Building MSA with {0} sequences.".format(len(related_samples)))

        # to just get the MSA
        msa_df = fn4c.msa(related_samples, output_format='json-records', what='N_or_M') 
        msa_df.set_index('guid',inplace = True)
        print(this_guid, msa_df.at[this_guid, 'p_value2'])


 