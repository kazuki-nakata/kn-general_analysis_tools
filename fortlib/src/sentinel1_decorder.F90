MODULE FDBAQDecoder
    IMPLICIT NONE
    INTEGER, PARAMETER :: BRC_BITS = 3
    INTEGER, PARAMETER :: THIDX_BITS = 8

    ! Huffman Trees for different BRC modes
    INTEGER, DIMENSION(0:4),PARAMETER :: MAX_HUFFMAN=(/3,4,6,9,15/)


CONTAINS

    ! メインのデコードサブルーチン
    SUBROUTINE DecodeFDBAQ(data,num_bytes,num_baq_blocks, num_quads, num_type, num_node, huffman_tree,& 
    &decoded_IE, decoded_IO, decoded_QE, decoded_QO,BRC,THIDX)
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: num_baq_blocks,num_quads,num_type,num_node,num_bytes
        ! CHARACTER(LEN=*), INTENT(IN) :: data
        INTEGER(1), DIMENSION(num_bytes), INTENT(IN) :: data
        INTEGER, DIMENSION(0:num_type-1, 0:num_node-1, 0:1), INTENT(IN) :: huffman_tree !tree_value:( +: node value -: minus decorded value )
        INTEGER, PARAMETER :: MAX_BLOCK_SIZE = 128
        INTEGER, DIMENSION(num_baq_blocks), INTENT(OUT) :: BRC, THIDX
        INTEGER, DIMENSION(MAX_BLOCK_SIZE, 2, num_baq_blocks), INTENT(OUT) :: decoded_IE, decoded_IO, decoded_QE, decoded_QO
        INTEGER :: block_index, bit_counter, byte_counter,i

        BRC(:)=0
        THIDX(:)=0
        decoded_IE(:,:,:)=0
        decoded_IO(:,:,:)=0
        decoded_QE(:,:,:)=0
        decoded_QO(:,:,:)=0

        ! バイト・ビットカウンタの初期化
        byte_counter = 1
        bit_counter = 0

        ! IE (Even In-phase) の処理
        ! ---------------------
        CALL ProcessChannelwithBRC(data, num_bytes, num_baq_blocks, num_quads, max_block_size, num_type,num_node, &
        &huffman_tree, decoded_IE, byte_counter, bit_counter, BRC)
    
        ! ! 16ビット境界への調整
        CALL AlignBit(byte_counter, bit_counter)

        ! ! ---------------------
        ! ! IO (Odd In-phase) の処理
        ! ! ---------------------
        CALL ProcessChannel(data, num_bytes, num_baq_blocks, num_quads, max_block_size, num_type,num_node, &
        &huffman_tree, decoded_IO, byte_counter, bit_counter, BRC)

        ! ! 16ビット境界への調整
        CALL AlignBit(byte_counter, bit_counter)

        ! ! ---------------------
        ! ! QE (Even Quadrature) の処理（THIDXを読む）
        ! ! ---------------------
        CALL ProcessChannelwithTHIDX(data, num_bytes, num_baq_blocks, num_quads, max_block_size, num_type, &
        &num_node,huffman_tree, decoded_QE, byte_counter, bit_counter, BRC, THIDX)

        ! ! 16ビット境界への調整
        CALL AlignBit(byte_counter, bit_counter)

        ! ! ---------------------
        ! ! QO (Odd Quadrature) の処理
        ! ! ---------------------
        CALL ProcessChannel(data, num_bytes, num_baq_blocks, num_quads, max_block_size, num_type,num_node, &
        &huffman_tree, decoded_QO, byte_counter, bit_counter, BRC)

    END SUBROUTINE DecodeFDBAQ

    SUBROUTINE ProcessChannelwithBRC(data, num_bytes, num_baq_blocks, num_quads, max_block_size, num_type, num_node, huffman_tree, &
    &decoded_data, byte_counter, bit_counter, BRC)
        INTEGER, INTENT(IN) :: num_bytes
        ! CHARACTER(LEN=*), INTENT(IN) :: data
        INTEGER(1), DIMENSION(num_bytes), INTENT(IN) :: data
        INTEGER, INTENT(IN) :: num_baq_blocks,max_block_size,num_quads,num_node,num_type
        INTEGER, DIMENSION(max_block_size, 2, num_baq_blocks), INTENT(OUT) :: decoded_data
        INTEGER, DIMENSION(num_baq_blocks), INTENT(OUT) :: BRC
        INTEGER, INTENT(INOUT) :: byte_counter, bit_counter
        INTEGER :: block_index, i, sign, current_node, decoded_value, outbit, block_size
        INTEGER, DIMENSION(0:num_type-1, 0:num_node-1, 0:1), INTENT(IN) :: huffman_tree !tree_value:( +: node value -: minus decorded value )
        decoded_data(:,:,:)=-999

        DO block_index = 1,num_baq_blocks
            CALL ReadBits(data, num_bytes, byte_counter, bit_counter, outbit, BRC_BITS)
            BRC(block_index) = outbit
            if(block_index.eq.num_baq_blocks) then
                block_size=num_quads-max_block_size*(num_baq_blocks-1)
            else
                block_size=128
            endif
            DO i = 1, block_size
                 CALL ReadBits(data, num_bytes, byte_counter, bit_counter, sign, 1)
                 CALL DecodeHuffman(data, num_bytes, byte_counter, bit_counter,num_node,&
                 &huffman_tree(BRC(block_index),:,:), decoded_value)
                 decoded_data(i, 1, block_index)= sign
                 decoded_data(i, 2, block_index)= decoded_value! * (2 * sign - 1)
            END DO
        END DO
    END SUBROUTINE ProcessChannelwithBRC

    SUBROUTINE ProcessChannel(data, num_bytes, num_baq_blocks, num_quads, max_block_size, num_type, num_node, huffman_tree, &
    &decoded_data, byte_counter, bit_counter, BRC)
        INTEGER, INTENT(IN) :: num_bytes
        ! CHARACTER(LEN=*), INTENT(IN) :: data
        INTEGER(1), DIMENSION(num_bytes), INTENT(IN) :: data
        INTEGER, INTENT(IN) :: num_baq_blocks,max_block_size,num_quads,num_node,num_type
        INTEGER, DIMENSION(max_block_size, 2, num_baq_blocks), INTENT(OUT) :: decoded_data
        INTEGER, DIMENSION(num_baq_blocks), INTENT(IN) :: BRC
        INTEGER, INTENT(INOUT) :: byte_counter, bit_counter
        INTEGER :: block_index, i, sign, current_node, decoded_value, outbit, block_size
        INTEGER, DIMENSION(0:num_type-1, 0:num_node-1, 0:1), INTENT(IN) :: huffman_tree !tree_value:( +: node value -: minus decorded value )
        decoded_data(:,:,:)=-999

        DO block_index = 1,num_baq_blocks
            if(block_index.eq.num_baq_blocks) then
                block_size=num_quads-max_block_size*(num_baq_blocks-1)
                ! block_size=max_block_size*num_baq_blocks-num_quads
            else
                block_size=128
            endif
            DO i = 1, block_size
                 CALL ReadBits(data, num_bytes, byte_counter, bit_counter, sign, 1)
                 CALL DecodeHuffman(data, num_bytes, byte_counter, bit_counter,num_node,&
                 &huffman_tree(BRC(block_index),:,:), decoded_value)
                 decoded_data(i, 1, block_index)= sign
                 decoded_data(i, 2, block_index)= decoded_value! * (2 * sign - 1)
            END DO
        END DO
    END SUBROUTINE ProcessChannel


    SUBROUTINE ProcessChannelwithTHIDX(data, num_bytes, num_baq_blocks, num_quads, max_block_size, num_type, num_node, &
    &huffman_tree, decoded_data, byte_counter, bit_counter, BRC, THIDX)
        INTEGER, INTENT(IN) :: num_bytes
        ! CHARACTER(LEN=*), INTENT(IN) :: data
        INTEGER(1), DIMENSION(num_bytes), INTENT(IN) :: data
        INTEGER, INTENT(IN) :: num_baq_blocks,max_block_size,num_quads,num_node,num_type
        INTEGER, DIMENSION(max_block_size, 2, num_baq_blocks), INTENT(OUT) :: decoded_data
        INTEGER, DIMENSION(num_baq_blocks), INTENT(IN) :: BRC
        INTEGER, DIMENSION(num_baq_blocks), INTENT(OUT) :: THIDX
        INTEGER, INTENT(INOUT) :: byte_counter, bit_counter
        INTEGER :: block_index, i, sign, current_node, decoded_value, outbit, block_size
        INTEGER, DIMENSION(0:num_type-1, 0:num_node-1, 0:1), INTENT(IN) :: huffman_tree !tree_value:( +: node value -: minus decorded value )
        decoded_data(:,:,:)=-999

        DO block_index = 1,num_baq_blocks
            CALL ReadBits(data, num_bytes, byte_counter, bit_counter, outbit, THIDX_BITS)
            THIDX(block_index) = outbit
            if(block_index.eq.num_baq_blocks) then
                block_size=num_quads-max_block_size*(num_baq_blocks-1)
            else
                block_size=128
            endif
            DO i = 1, block_size
                 CALL ReadBits(data, num_bytes, byte_counter, bit_counter, sign, 1)
                 CALL DecodeHuffman(data, num_bytes, byte_counter, bit_counter,num_node,&
                 &huffman_tree(BRC(block_index),:,:), decoded_value)
                 decoded_data(i, 1, block_index)= sign
                 decoded_data(i, 2, block_index)= decoded_value! * (2 * sign - 1)
            END DO
        END DO
    END SUBROUTINE ProcessChannelwithTHIDX


    ! THIDX を読み取るサブルーチン
    SUBROUTINE ReadTHIDX(data, num_bytes, num_baq_blocks, THIDX, byte_counter, bit_counter)
        INTEGER, INTENT(IN) :: num_bytes
        INTEGER(1), DIMENSION(num_bytes), INTENT(IN) :: data
        INTEGER, INTENT(IN) :: num_baq_blocks
        INTEGER, DIMENSION(:), INTENT(OUT) :: THIDX
        INTEGER, INTENT(INOUT) :: byte_counter, bit_counter
        INTEGER :: block_index,outbit

        DO block_index = 1, num_baq_blocks
            CALL ReadBits(data, num_bytes, byte_counter, bit_counter, outbit, THIDX_BITS)
            THIDX(block_index) = outbit
        END DO
    END SUBROUTINE ReadTHIDX

    ! 1ビットずつ読み取る
    SUBROUTINE ReadBits(data, num_bytes, byte_counter, bit_counter, value, num_bits)
        INTEGER, INTENT(IN) :: num_bytes
        INTEGER(1), DIMENSION(num_bytes), INTENT(IN) :: data
        ! CHARACTER(LEN=*), INTENT(IN) :: bit_data
        INTEGER, INTENT(INOUT) :: byte_counter, bit_counter
        INTEGER, INTENT(IN) :: num_bits
        INTEGER, INTENT(OUT) :: value
        INTEGER(4) :: ascii
        INTEGER(4) :: i
        value = 0

        DO i = 1, num_bits
            ascii= data(byte_counter)
            value = IBITS(ascii, 7-bit_counter, 1) + 2 * value
            bit_counter = bit_counter + 1
            IF (bit_counter == 8) THEN
                bit_counter = 0
                byte_counter = byte_counter + 1
            END IF
        END DO
    END SUBROUTINE ReadBits

    ! 16ビット境界への調整
    SUBROUTINE AlignBit(byte_counter, bit_counter)
        INTEGER, INTENT(INOUT) :: byte_counter, bit_counter
        IF (bit_counter /= 0) THEN
            bit_counter = 0
            byte_counter = byte_counter + 1
        END IF
        byte_counter = int(real(byte_counter) / 2) * 2 + 1
    END SUBROUTINE AlignBit

   ! Huffman 復号
    SUBROUTINE DecodeHuffman(endata, num_bytes, byte_counter, bit_counter, num_node, huffman_tree, decoded_value) 
        ! CHARACTER(LEN=*), INTENT(IN) :: endata
        INTEGER, INTENT(IN) :: num_bytes
        INTEGER(1), DIMENSION(num_bytes), INTENT(IN) :: endata
        INTEGER, INTENT(INOUT) :: byte_counter, bit_counter
        INTEGER, INTENT(OUT) :: decoded_value
        INTEGER, INTENT(IN) :: num_node
        INTEGER :: current_node,outbit,i,max_huffman0,value
        INTEGER, DIMENSION(0:num_node-1,0:1), INTENT(IN) :: huffman_tree
        current_node=0
        value=1
        DO WHILE (value>0.5)
            CALL ReadBits(endata, num_bytes, byte_counter, bit_counter, outbit, 1)
            current_node=huffman_tree(current_node,outbit)
            value=huffman_tree(current_node,outbit)
        END DO
        decoded_value = -value
    END SUBROUTINE DecodeHuffman

    SUBROUTINE getSentinel1RawMap(data, brc, thidx, max_block_size, num_baq_blocks,b0,b1,b2,b3,b4,&
    &nrl_b0,nrl_b1,nrl_b2,nrl_b3,nrl_b4,sf,output)
        INTEGER, INTENT(IN) :: max_block_size,num_baq_blocks
        INTEGER(4), DIMENSION(max_block_size, 0:1, num_baq_blocks), INTENT(IN) :: data
        INTEGER(4), DIMENSION(num_baq_blocks), INTENT(IN) :: brc, thidx
        REAL(4), DIMENSION(max_block_size,num_baq_blocks), INTENT(OUT) :: output
        INTEGER(4) :: i,j
        INTEGER(4), DIMENSION(0:4) :: thidx_thres,data_thres
        REAL(4) :: b0(0:3),b1(0:3),b2(0:5),b3(0:6),b4(0:8)
        REAL(4) :: nrl_b0(0:3),nrl_b1(0:4),nrl_b2(0:6),nrl_b3(0:9),nrl_b4(0:15),sf(0:255)
        REAL(4),DIMENSION(max_block_size) :: coef1,coef2
        
        thidx_thres=(/3,3,5,6,8/)
        data_thres=(/3,4,6,9,15/)
        output(:,:)=0

        do j = 1, num_baq_blocks
            if(brc(j).eq.0) then
                if(thidx(j).le.thidx_thres(0)) then
                    coef2(:)=1
                    do i =1,max_block_size
                        if(data(i,1,j).lt.data_thres(0)) then
                            coef1(i)=data(i,1,j)
                        else
                            coef1(i)=b0(thidx(j))
                        endif
                    enddo                            
                else
                    coef2(:)=sf(thidx(j))
                    do i =1,max_block_size
                        coef1(i)=nrl_b0(data(i,1,j))
                    enddo
                endif
            elseif(brc(j).eq.1) then
                if(thidx(j).le.thidx_thres(1)) then
                    coef2(:)=1
                    do i =1,max_block_size
                        if(data(i,1,j).lt.data_thres(1)) then
                            coef1(i)=data(i,1,j)
                        else
                            coef1(i)=b1(thidx(j))
                        endif
                    enddo                            
                else
                    coef2(:)=sf(thidx(j))
                    do i =1,max_block_size
                        coef1(i)=nrl_b1(data(i,1,j))
                    enddo
                endif
            elseif(brc(j).eq.2) then
                if(thidx(j).le.thidx_thres(2)) then
                    coef2(:)=1
                    do i =1,max_block_size
                        if(data(i,1,j).lt.data_thres(2)) then
                            coef1(i)=data(i,1,j)
                        else
                            coef1(i)=b2(thidx(j))
                        endif
                    enddo                            
                else
                    coef2(:)=sf(thidx(j))
                    do i =1,max_block_size
                        coef1(i)=nrl_b2(data(i,1,j))
                    enddo
                endif
            elseif(brc(j).eq.3) then
                if(thidx(j).le.thidx_thres(3)) then
                    coef2(:)=1
                    do i =1,max_block_size
                        if(data(i,1,j).lt.data_thres(3)) then
                            coef1(i)=data(i,1,j)
                        else
                            coef1(i)=b3(thidx(j))
                        endif
                    enddo                            
                else
                    coef2(:)=sf(thidx(j))
                    do i =1,max_block_size
                        coef1(i)=nrl_b3(data(i,1,j))
                    enddo
                endif
            elseif(brc(j).eq.4) then
                if(thidx(j).le.thidx_thres(4)) then
                    coef2(:)=1
                    do i =1,max_block_size
                        if(data(i,1,j).lt.data_thres(4)) then
                            coef1(i)=data(i,1,j)
                        else
                            coef1(i)=b4(thidx(j))
                        endif
                    enddo                            
                else
                    coef2(:)=sf(thidx(j))
                    do i =1,max_block_size
                        coef1(i)=nrl_b4(data(i,1,j))
                    enddo
                endif
            endif
            do i = 1, max_block_size
                output(i,j)=((-1)**data(i,0,j))*coef1(i)*coef2(i)
            enddo
        enddo
    END SUBROUTINE getSentinel1RawMap

    SUBROUTINE AllProc(data,num_bytes0,num_line,num_baq_blocks, num_quads, num_bytes,num_type,&
    & num_node, huffman_tree, param_dir,output)
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: num_baq_blocks,num_quads,num_type,num_node,num_bytes0,num_line
        INTEGER, DIMENSION(num_line) :: num_bytes
        CHARACTER(LEN=*), INTENT(IN) :: param_dir
        INTEGER(1), DIMENSION(num_bytes0,num_line), INTENT(IN) :: data
        INTEGER, DIMENSION(0:num_type-1, 0:num_node-1, 0:1), INTENT(IN) :: huffman_tree
        COMPLEX(4), DIMENSION(num_quads*2,num_line), INTENT(OUT) :: output
        INTEGER, PARAMETER :: MAX_BLOCK_SIZE = 128
        INTEGER, DIMENSION(num_baq_blocks) :: BRC, THIDX
        INTEGER, DIMENSION(MAX_BLOCK_SIZE, 0:1, num_baq_blocks) :: decoded_IE, decoded_IO, decoded_QE, decoded_QO
        REAL(4), DIMENSION(MAX_BLOCK_SIZE, num_baq_blocks) :: IE, IO, QE, QO
        INTEGER(1), DIMENSION(MAX_BLOCK_SIZE, num_baq_blocks) :: mask
        INTEGER :: i,j,i2,j2
        REAL(4) :: b0(0:3),b1(0:3),b2(0:5),b3(0:6),b4(0:8)
        REAL(4) :: nrl_b0(0:3),nrl_b1(0:4),nrl_b2(0:6),nrl_b3(0:9),nrl_b4(0:15),sf(0:255)

        open(21,file=TRIM(param_dir)//'lookup_b0.dat', status='old',form='unformatted',access='direct',recl=4*4)
        read(21,rec=1) b0
        close(21)
        open(21,file=TRIM(param_dir)//'lookup_b1.dat', status='old',form='unformatted',access='direct',recl=4*4)
        read(21,rec=1) b1
        close(21)
        open(21,file=TRIM(param_dir)//'lookup_b2.dat', status='old',form='unformatted',access='direct',recl=4*6)
        read(21,rec=1) b2
        close(21)
        open(21,file=TRIM(param_dir)//'lookup_b3.dat', status='old',form='unformatted',access='direct',recl=4*7)
        read(21,rec=1) b3
        close(21)
        open(21,file=TRIM(param_dir)//'lookup_b4.dat', status='old',form='unformatted',access='direct',recl=4*9)
        read(21,rec=1) b4
        close(21)
        open(21,file=TRIM(param_dir)//'lookup_nrl_b0.dat', status='old',form='unformatted',access='direct',recl=4*4)
        read(21,rec=1) nrl_b0
        close(21)
        open(21,file=TRIM(param_dir)//'lookup_nrl_b1.dat', status='old',form='unformatted',access='direct',recl=4*5)
        read(21,rec=1) nrl_b1
        close(21)
        open(21,file=TRIM(param_dir)//'lookup_nrl_b2.dat', status='old',form='unformatted',access='direct',recl=4*7)
        read(21,rec=1) nrl_b2
        close(21)
        open(21,file=TRIM(param_dir)//'lookup_nrl_b3.dat', status='old',form='unformatted',access='direct',recl=4*10)
        read(21,rec=1) nrl_b3
        close(21)
        open(21,file=TRIM(param_dir)//'lookup_nrl_b4.dat', status='old',form='unformatted',access='direct',recl=4*16)
        read(21,rec=1) nrl_b4
        close(21)
        open(21,file=TRIM(param_dir)//'lookup_sf.dat', status='old',form='unformatted',access='direct',recl=4*256)
        read(21,rec=1) sf
        close(21)
        output(:,:)=0.
        do j = 1, num_line
            CALL DecodeFDBAQ(data(1:num_bytes(j),j),num_bytes(j),num_baq_blocks, num_quads, num_type, num_node, huffman_tree,& 
                &decoded_IE, decoded_IO, decoded_QE, decoded_QO,BRC,THIDX)
            
            where(decoded_IE(:,0,:).eq.-999)
                mask=0
            else where
                mask=1
            end where
            where(decoded_IE.eq.-999) decoded_IE=0
            where(decoded_IO.eq.-999) decoded_IO=0
            where(decoded_QE.eq.-999) decoded_QE=0
            where(decoded_QO.eq.-999) decoded_QO=0
            CALL getSentinel1RawMap(decoded_IE, brc, thidx, max_block_size, num_baq_blocks,b0,b1,b2,b3,b4,&
    &nrl_b0,nrl_b1,nrl_b2,nrl_b3,nrl_b4,sf,IE)
            CALL getSentinel1RawMap(decoded_IO, brc, thidx, max_block_size, num_baq_blocks,b0,b1,b2,b3,b4,&
    &nrl_b0,nrl_b1,nrl_b2,nrl_b3,nrl_b4,sf,IO)
            CALL getSentinel1RawMap(decoded_QE, brc, thidx, max_block_size, num_baq_blocks,b0,b1,b2,b3,b4,&
    &nrl_b0,nrl_b1,nrl_b2,nrl_b3,nrl_b4,sf,QE)
            CALL getSentinel1RawMap(decoded_IE, brc, thidx, max_block_size, num_baq_blocks,b0,b1,b2,b3,b4,&
    &nrl_b0,nrl_b1,nrl_b2,nrl_b3,nrl_b4,sf,QO)
            i=1
            do j2=1,num_baq_blocks
                do i2=1,max_block_size
                    if(mask(i2,j2).eq.1) then
                        output(i,j)= cmplx(IE(i2,j2),QE(i2,j2))
                        i=i+1
                        output(i,j)= cmplx(IO(i2,j2),QO(i2,j2))
                        i=i+1
                    endif
                enddo
            enddo
        enddo

    ENDSUBROUTINE AllProc


END MODULE FDBAQDecoder