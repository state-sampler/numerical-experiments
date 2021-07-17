function random_int(start, stop) result(i)
    implicit none
    integer, intent(in) :: start, stop
    real*8  :: u
    integer :: i

    call random_number(u)
    i = start + FLOOR((stop + 1 - start)*u)
end function random_int

subroutine gibbs_sampler(set, partial_order, chain, chain_index, Ei, card_diff, index_size, size)
  ! Basic implementation of a Gibbs sampler.  Given a set in E_i, the sampler returns another set in E_i by block sampling on the chains extracted from the partial order.
  implicit none
  integer*8, dimension(size), intent(inout) :: set ! Starting set of items that is updated with the latest sample set
  integer*8, dimension(size,size), intent(in) :: partial_order ! Partial_order(i, j) == 1 iff i <= j in the partial order
  integer*8, dimension(size), intent(in) :: chain ! Array containing the concatenated chains
  integer*8, dimension(index_size), intent(in) :: chain_index ! Array containing the starting index for each chain in the chain array
  integer*8, intent(in) :: size, index_size, Ei ! Total number of items, total number of chains, current value of E_i
  integer*8, dimension(1), intent(inout) :: card_diff ! Distance from the current set to the smallest knowledge state containing it

  ! Miscellaneous variables used when running the sampler
  integer*8, dimension(size) :: closed_set, curr_set, ind_list, chain_diff
  integer :: i, j, k, n, added_item, count, chain_size, start_ind, end_ind, ind
  integer, external :: random_int

  ! Main loop that (block) samples from within each of the chains
  do n=1,index_size
     ! Find starting and ending indices for current chain
     start_ind=chain_index(n)
     if (n==index_size) then
        end_ind=size
     else
        end_ind=chain_index(n+1)-1
     endif
     chain_size=1+end_ind-start_ind

     ! Initialize array for computing d_F(X) values (i.e., for each set X, the distance from the smallest knowledge state containing X)
     do j=1,chain_size+1
        chain_diff(j)=0
     enddo
     ! Remove the chain from the current set
     do j=start_ind,end_ind
        set(chain(j))=0
     enddo
     ! Initialize arrays for subsequent computations
     do j=1,size
        closed_set(j)=set(j)
        curr_set(j)=set(j)
     enddo

     ! Find the lower closure of the current set (i.e., for each item in the current set, add all the items that come before it in the partial order)
     do j=1,size
        if (curr_set(j)==1) then
           do i=1,size
              if (partial_order(i,j)==1) then
                 closed_set(i)=1
              endif
           enddo
        endif
     enddo

     ! Compute distance from the current set to its lower closure
     do i=1,size
        chain_diff(1)=chain_diff(1)+(closed_set(i)-curr_set(i))
     enddo
    ! Iteratively add items from chain to current set and compute distance to lower closure
     do k=1,chain_size
        added_item=chain(k-1+start_ind)
        curr_set(added_item)=1
        do i=1,size
           if (partial_order(i,added_item)==1) then
              closed_set(i)=1
           endif
        enddo
        do i=1,size
           chain_diff(k+1)=chain_diff(k+1)+(closed_set(i)-curr_set(i))
        enddo
     enddo
     ! Find sets contained in E_i
     count=0
     do k=1,chain_size+1
        if(chain_diff(k)<=Ei) then
           count=count+1
           ind_list(count)=k
        endif
     enddo
     ! Randomly choose one of the sets in E_i 
     ind=ind_list(random_int(1,count))
     if(ind>=2) then
        do i=start_ind,start_ind+ind-2
           set(chain(i))=1
        enddo
     endif
     if (n==index_size) then
        card_diff(1)=chain_diff(ind)
     endif
  enddo
end subroutine gibbs_sampler
