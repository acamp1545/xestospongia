UPDATE data
SET Kingdom=0
WHERE Kingdom IS NULL;

UPDATE data
SET Phylum=0
WHERE Phylum IS NULL;

UPDATE data
SET Class=0
WHERE Class IS NULL;

UPDATE data
SET `Order`=0
WHERE `Order` IS NULL;

UPDATE data
SET Family=0
WHERE Family IS NULL;

UPDATE data
SET Genus=0
WHERE Genus IS NULL;

UPDATE data
SET Species=0
WHERE Species IS NULL;