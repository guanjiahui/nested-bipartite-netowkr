% Several examples, using the Matlab wrapper.

examples_to_run = [1.1,1.2,2.1,3.1,3.2];
number_of_samples = 10;

for example = examples_to_run
    switch example
        % ==================================================================
        % Ecology co-occurrence matrices
        case 1.1
        title = 'Darwin''s finches';
        p = [14 13 14 10 12 2 10 1 10 11 6 2 17];
        q = [4 4 11 10 10 8 9 10 8 9 3 10 4 7 9 3 3];
        matrix_type = 0;
        true_number = '67149106137567626';

        case 1.2
        title = 'Manly''s lizards';
        p = [7 23 6 6 2 18 8 8 1 22 2 2 9 2 1 18 9 3 3 1];
        q = [13 4 9 4 3 2 3 4 4 5 2 5 2 10 10 10 7 6 6 3 3 11 8 11 6];
        matrix_type = 0;
        true_number = '55838420515731001979319625577023858901579264';
        
        case 1.3
        title = 'Patterson & Atmar''s montane mammals';
        p = [26 26 25 22 22 18 12 12 12 11 10 10 8 8 8 7 6 6 5 5 4 4 3 3 1 1];
        q = [26 24 23 21 19 13 13 12 11 10 10 9 9 7 7 7 7 7 7 6 6 5 5 4 3 2 1 1];
        matrix_type = 0;
        true_number = '2663296694330271332856672902543209853700';
        
        % ==================================================================
        % Social networks
        case 2.1
        title = 'Galaskiewicz''s CEOs and clubs';
        p = [3 11 22 12 3 4 4 4 6 3 4 5 5 3 9];
        q = [3 3 2 3 3 3 4 3 4 2 3 2 4 7 5 5 6 5 5 5 3 3 4 5 3 3]; 
        matrix_type = 0;
        true_number = '25533540876226059861329182058955286218365274646344655';
        
        case 2.2
        title = 'Ragusan marriages';
        p = [3 0 3 0 19 0 4 7 5 4 12 7 4 8 0 4 1 1 2 5 0 14 2 8]; 
        q = [0 1 4 2 20 2 4 4 5 4 15 9 3 6 1 2 0 0 2 4 1 21 0 3];
        matrix_type = 1;
        true_number = '949599133340064609956529916243690765923314405123631501076903661';
        
        % ==================================================================
        % Contingency tables
        case 3.1
        title = 'Galton''s couples';
        p = [51 104 50];
        q = [46 99 60];
        matrix_type = 1;
        true_number = '1268792';

        case 3.2
        title = 'Diaconis & Gangolli, Example (1)';
        p = [65 25 45];
        q = [10 62 13 11 39];
        matrix_type = 1;
        true_number = '239382173';

        case 3.3
        title = 'Diaconis & Gangolli, Example (2)';
        p = [108  286  71  127];
        q = [220  215  93  64];
        matrix_type = 1;
        true_number = '1225914276768514';

        % ==================================================================
        % Large example
        case 4.1
        title = 'Large example (binary)';
        v = @(k) ones(1,k);
        p = [5*v(20),4*v(20),3*v(20),2*v(20),1*v(20)];
        q = p;
        matrix_type = 0;
        true_number = '235147658512972346136720844732752867234190813428532503765177603865861155465020348800377149750640574119630917726455436243685355314182055905496772141379315815531771443489007534912882316530744246557700549985656553879075329865611504724448050744461432151755297042930358350343972327816377263569621228846426652684966621966838583286516872063741754047640401277490098051066579339199616742757648874362368765426913127210701207633920000000000000';
        
        case 4.2
        title = 'Large example (nonnegative integer)';
        v = @(k) ones(1,k);
        p = [5*v(20),4*v(20),3*v(20),2*v(20),1*v(20)];
        q = p;
        matrix_type = 1;
        true_number = '295805666559153333218627349404825478412279341291256310789770124587755158139247610594820909327345130630237874869835355865067328480221768746326533318518034068762147953091161124607645707089033147454208952389634816899177414716448735101144730357394276262070385948234353936879761227595488660906921645496874333081165039663928985440675621736614359726081137694371480321973707301121622703911176113569810500814339034022648523951964160000000000000';
    end
    
    % Summary
    fprintf('\n============== Example %.1f ==============\n',example);
    fprintf([title '\n']);
    fprintf('p = row sums = [ '); fprintf('%d ',p); fprintf(']\n');
    fprintf('q = col sums = [ '); fprintf('%d ',q); fprintf(']\n');
    
    % Count the number of matrices with row sums p and column sums q
    fprintf('Counting...\n');
    tic
    number = count(p,q,matrix_type);
    toc
    % (or, to specify the filenames to be used:)
    % number = count(p,q,matrix_type,'input.txt','table.bin')

    % Verify that we got the correct answer
    fprintf('True =    %s\n',true_number);
    fprintf('Counted = %s\n',number);
    assert(strcmp(number,true_number));

    % Draw some random samples from the uniform distribution on this set of matrices
    fprintf('Sampling %d matrices...\n',number_of_samples);
    X = sample(p,q,number_of_samples);
    % (or, to specify the filenames to be used:)
    % X = sample(p,q,number_of_samples,'input.txt','table.bin','output.txt');
    
    % Display the 1st sample (use fprintf to avoid line wrapping)
    fprintf('1st sample =\n');
    for i = 1:length(p); fprintf('%d ',X(i,:,1)); fprintf('\n'); end;

    % Verify that the samples have the correct row and column sums
    for i = 1:number_of_samples
        assert(all(sum(X(:,:,i),1)==q));
        assert(all(sum(X(:,:,i),2)==p'));
    end
end



