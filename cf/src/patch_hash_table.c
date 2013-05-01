#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include "util.h"
#include "patch_hash_table.h"

TableKeyType _hash(KeyType originKey, size_t tableSize) {
#ifdef DEBUG
	fprintf(stderr, "hash: tableSize: %d\n", tableSize);
	fprintf(stderr, "hash: mult: %d\n", mult);
	fprintf(stderr, "hash: %d\n", hash);
#endif
	int mult = (originKey.patchID * originKey.hillID * originKey.zoneID);
	TableKeyType hash = (TableKeyType) (mult % tableSize);
	return hash;
}

bool _keysAreEqual(KeyType key1, KeyType key2) {
	if (key1.patchID != key2.patchID) return false;
	if (key1.hillID != key2.hillID) return false;
	if (key1.zoneID != key2.zoneID) return false;

	return true;
}

bool _keyIsEmpty(KeyType key) {
	if ((PATCH_HASH_TABLE_EMPTY == key.patchID)
			&& (PATCH_HASH_TABLE_EMPTY == key.hillID)
			&& (PATCH_HASH_TABLE_EMPTY == key.zoneID)) {
		return true;
	}

	return false;
}

Table *allocatePatchHashTable(size_t tableSize) {
	Table* table = (Table *) malloc(sizeof(Table));
	table->tableSize = tableSize;
	table->numEntries = 0;
	table->entries = (TableEntry *) calloc(tableSize, sizeof(TableEntry));
	if (NULL == table->entries) {
		free(table);
		table = NULL;
	}
	// Initialize keys to empty values
	TableEntry *entry;
	size_t size = sizeof(TableEntry);
	for (size_t i = 0; i < tableSize; i++) {
		entry = (TableEntry *) (table->entries + (i * size));
		entry->originKey.patchID = PATCH_HASH_TABLE_EMPTY;
		entry->originKey.hillID = PATCH_HASH_TABLE_EMPTY;
		entry->originKey.zoneID = PATCH_HASH_TABLE_EMPTY;
		entry->value = PATCH_HASH_TABLE_EMPTY;
		entry->next = NULL;
	}
	return table;
}

void freePatchHashTable(Table *table) {
	assert(table != NULL);
	// Free buckets
	TableEntry *entry;
	TableEntry *tmpEntry;
	size_t size = sizeof(TableEntry);
	for (size_t i = 0; i < table->tableSize; i++) {
		entry = (TableEntry *) (table->entries + (i * size));
		tmpEntry = entry->next;
		while (tmpEntry != NULL ) {
			entry = tmpEntry;
			tmpEntry = entry->next;
			free(entry);
		}
	}
	free(table->entries);
	free(table);
}

void patchHashTableInsert(Table *table, KeyType key, ValueType value) {
	assert(table != NULL);
	TableKeyType hash = _hash(key, table->tableSize);
#ifdef DEBUG
	fprintf(stderr, "Inserting key: %d, %d, %d; value: %d; hash: %d\n",
				key.patchID, key.hillID, key.zoneID, value, hash);
#endif
	TableEntry *entry = (TableEntry *) ( table->entries + ( hash * sizeof(TableEntry) ) );
	if (_keyIsEmpty(entry->originKey)) {
		// Insert first key
#ifdef DEBUG
		fprintf(stderr, "key is empty\n");
#endif
		entry->originKey = key;
		entry->value = value;
		entry->next = NULL;
		table->numEntries++;
	} else {
		while (entry) {
			if (_keysAreEqual(entry->originKey, key)) {
				// Key is already present, overwrite value and return
#ifdef DEBUG
				fprintf(stderr, "keys are equal\n");
#endif
				entry->value = value;
				return;
			}
			if (NULL == entry->next) {
				// We've reached the end of the list, add a new entry and return
#ifdef DEBUG
				fprintf(stderr, "end of list; next: %p\n", entry->next);
#endif
				entry->next = (TableEntry *)malloc( sizeof(TableEntry) );
#ifdef DEBUG
				fprintf(stderr, "next: %p\n", entry->next);
#endif
				assert( entry->next );
				entry->next->originKey = key;
				entry->next->value = value;
				entry->next->next = NULL;
				table->numEntries++;
				return;
			}
			entry = entry->next;
		}
	}
}

ValueType patchHashTableGet(Table *table, KeyType key) {
	assert(table != NULL);
	TableKeyType hash = _hash(key, table->tableSize);
	TableEntry *entry = (TableEntry *) (table->entries
			+ (hash * sizeof(TableEntry)));
	// Search for key
	while (entry) {
		if (_keysAreEqual(entry->originKey, key)) {
#ifdef DEBUG
			fprintf(stderr, "found key: %d, %d, %d, returning value: %d\n",
					key.patchID, key.hillID, key.zoneID, entry->value);
#endif
			return entry->value;
		}
		entry = entry->next;
	}
	return PATCH_HASH_TABLE_EMPTY;
}

void printPatchHashTable(Table *table) {
	TableEntry *entry;
	TableEntry *tmpEntry;
	size_t size = sizeof(TableEntry);
	for (size_t i = 0; i < table->tableSize; i++) {
		entry = (TableEntry *) (table->entries + (i * size));
		fprintf(stderr, "Entry %d: originKey: %d, %d, %d; value: %d; next: %p\n", i,
				entry->originKey.patchID, entry->originKey.hillID,
				entry->originKey.zoneID, entry->value, entry->next);
		tmpEntry = entry->next;
		while (tmpEntry != NULL ) {
			entry = tmpEntry;

			fprintf(stderr, "\tEntry %d: originKey: %d, %d, %d; value: %d; next: %p\n", i,
					entry->originKey.patchID, entry->originKey.hillID,
					entry->originKey.zoneID, entry->value, entry->next);

			tmpEntry = entry->next;

		}
	}
}
