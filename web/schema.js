const schema = {
  $schema: "http://json-schema.org/draft-07/schema#",
  additionalProperties: false,
  definitions: {
    AtomJ1L2ConfigImpl: {
      additionalProperties: false,
      properties: {
        core: {
          $ref: "#/definitions/AtomLSImpl",
        },
        excited: {
          $ref:
            "#/definitions/ConfigTerm<CouplingScheme.LS,LSTermImpl,ShellEntry[]>",
        },
      },
      required: ["core", "excited"],
      type: "object",
    },
    AtomLSImpl: {
      additionalProperties: false,
      properties: {
        config: {
          items: {
            $ref: "#/definitions/ShellEntry",
          },
          type: "array",
        },
        scheme: {
          $ref: "#/definitions/CouplingScheme.LS",
        },
        term: {
          additionalProperties: false,
          properties: {
            J: {
              type: "number",
            },
            L: {
              type: "number",
            },
            P: {
              type: "number",
            },
            S: {
              type: "number",
            },
          },
          required: ["J", "L", "P", "S"],
          type: "object",
        },
      },
      required: ["config", "scheme", "term"],
      type: "object",
    },
    "ConfigTerm<CouplingScheme.LS,LSTermImpl,ShellEntry[]>": {
      additionalProperties: false,
      properties: {
        config: {
          items: {
            $ref: "#/definitions/ShellEntry",
          },
          type: "array",
        },
        scheme: {
          $ref: "#/definitions/CouplingScheme.LS",
        },
        term: {
          $ref: "#/definitions/LSTermImpl",
        },
      },
      required: ["config", "scheme", "term"],
      type: "object",
    },
    "CouplingScheme.J1L2": {
      enum: [1],
      type: "number",
    },
    "CouplingScheme.LS": {
      enum: [0],
      type: "number",
    },
    CrossSection: {
      additionalProperties: false,
      properties: {
        data: {
          items: {
            items: {
              type: "number",
            },
            type: "array",
          },
          type: "array",
        },
        labels: {
          items: {
            type: "string",
          },
          type: "array",
        },
        parameters: {
          items: {
            type: "string",
          },
          type: "array",
        },
        reaction: {
          $ref: "#/definitions/Reaction",
        },
        reference: {
          anyOf: [
            {
              items: {
                type: "string",
              },
              type: "array",
            },
            {
              items: {
                $ref: "#/definitions/Reference",
              },
              type: "array",
            },
          ],
        },
        threshold: {
          type: "number",
        },
        units: {
          items: {
            type: "string",
          },
          type: "array",
        },
      },
      required: ["data", "labels", "reaction", "threshold", "units"],
      type: "object",
    },
    LSTermImpl: {
      additionalProperties: false,
      properties: {
        L: {
          type: "number",
        },
        P: {
          type: "number",
        },
        S: {
          type: "number",
        },
      },
      required: ["L", "P", "S"],
      type: "object",
    },
    Reaction: {
      additionalProperties: false,
      properties: {
        lhs: {
          items: {
            $ref: "#/definitions/ReactionEntry",
          },
          type: "array",
        },
        reversible: {
          type: "boolean",
        },
        rhs: {
          items: {
            $ref: "#/definitions/ReactionEntry",
          },
          type: "array",
        },
        type_tags: {
          items: {
            enum: [
              "Effective",
              "Elastic",
              "Electronic",
              "Ionization",
              "Rotational",
              "Vibrational",
            ],
            type: "string",
          },
          type: "array",
        },
      },
      required: ["lhs", "reversible", "rhs", "type_tags"],
      type: "object",
    },
    ReactionEntry: {
      additionalProperties: false,
      properties: {
        count: {
          type: "number",
        },
        state: {
          type: "string",
        },
      },
      required: ["count", "state"],
      type: "object",
    },
    Reference: {
      additionalProperties: false,
      properties: {
        Email: {
          items: {
            type: "string",
          },
          type: "array",
        },
        address: {
          type: "string",
        },
        alt_type: {
          type: "string",
        },
        annote: {
          type: "string",
        },
        author: {
          items: {
            type: "string",
          },
          type: "array",
        },
        booktitle: {
          type: "string",
        },
        chapter: {
          type: "number",
        },
        crossref: {
          type: "string",
        },
        doi: {
          type: "string",
        },
        edition: {
          type: "string",
        },
        howpublished: {
          type: "string",
        },
        id: {
          type: "string",
        },
        institution: {
          type: "string",
        },
        journal: {
          type: "string",
        },
        key: {
          type: "string",
        },
        month: {
          type: "string",
        },
        note: {
          type: "string",
        },
        number: {
          type: "number",
        },
        organization: {
          type: "string",
        },
        pages: {
          type: "string",
        },
        publisher: {
          type: "string",
        },
        school: {
          type: "string",
        },
        series: {
          type: "string",
        },
        title: {
          type: "string",
        },
        type: {
          $ref: "#/definitions/ReferenceType",
        },
        volume: {
          type: "number",
        },
        year: {
          type: "number",
        },
      },
      required: ["author", "id", "title", "type", "year"],
      type: "object",
    },
    ReferenceType: {
      enum: [
        "article",
        "book",
        "booklet",
        "conference",
        "inbook",
        "incollection",
        "inproceedings",
        "manual",
        "mastersthesis",
        "misc",
        "phdthesis",
        "proceedings",
        "techreport",
        "unpublished",
      ],
      type: "string",
    },
    ShellEntry: {
      additionalProperties: false,
      properties: {
        l: {
          type: "number",
        },
        n: {
          type: "number",
        },
        occupance: {
          type: "number",
        },
      },
      required: ["l", "n", "occupance"],
      type: "object",
    },
  },
  properties: {
    processes: {
      items: {
        $ref: "#/definitions/CrossSection",
      },
      type: "array",
    },
    states: {
      items: {
        additionalProperties: false,
        properties: {
          charge: {
            type: "number",
          },
          electronic: {
            anyOf: [
              {
                items: {
                  anyOf: [
                    {
                      additionalProperties: false,
                      properties: {
                        config: {
                          items: {
                            $ref: "#/definitions/ShellEntry",
                          },
                          type: "array",
                        },
                        scheme: {
                          $ref: "#/definitions/CouplingScheme.LS",
                        },
                        summary: {
                          type: "string",
                        },
                        term: {
                          additionalProperties: false,
                          properties: {
                            J: {
                              type: "number",
                            },
                            L: {
                              type: "number",
                            },
                            P: {
                              type: "number",
                            },
                            S: {
                              type: "number",
                            },
                          },
                          required: ["J", "L", "P", "S"],
                          type: "object",
                        },
                      },
                      required: ["config", "scheme", "summary", "term"],
                      type: "object",
                    },
                    {
                      additionalProperties: false,
                      properties: {
                        e: {
                          type: "string",
                        },
                      },
                      required: ["e"],
                      type: "object",
                    },
                  ],
                },
                type: "array",
              },
              {
                items: {
                  anyOf: [
                    {
                      additionalProperties: false,
                      properties: {
                        config: {
                          $ref: "#/definitions/AtomJ1L2ConfigImpl",
                        },
                        scheme: {
                          $ref: "#/definitions/CouplingScheme.J1L2",
                        },
                        summary: {
                          type: "string",
                        },
                        term: {
                          additionalProperties: false,
                          properties: {
                            J: {
                              type: "number",
                            },
                            K: {
                              type: "number",
                            },
                            P: {
                              type: "number",
                            },
                            S: {
                              type: "number",
                            },
                          },
                          required: ["J", "K", "P", "S"],
                          type: "object",
                        },
                      },
                      required: ["config", "scheme", "summary", "term"],
                      type: "object",
                    },
                    {
                      additionalProperties: false,
                      properties: {
                        e: {
                          type: "string",
                        },
                      },
                      required: ["e"],
                      type: "object",
                    },
                  ],
                },
                type: "array",
              },
              {
                items: {
                  anyOf: [
                    {
                      additionalProperties: false,
                      properties: {
                        e: {
                          type: "string",
                        },
                        summary: {
                          type: "string",
                        },
                      },
                      required: ["e", "summary"],
                      type: "object",
                    },
                    {
                      additionalProperties: false,
                      properties: {
                        Lambda: {
                          type: "number",
                        },
                        S: {
                          type: "number",
                        },
                        e: {
                          type: "string",
                        },
                        parity: {
                          enum: ["g", "u"],
                          type: "string",
                        },
                        reflection: {
                          enum: ["+", "-"],
                          type: "string",
                        },
                        summary: {
                          type: "string",
                        },
                        vibrational: {
                          items: [
                            {
                              anyOf: [
                                {
                                  additionalProperties: false,
                                  properties: {
                                    summary: {
                                      type: "string",
                                    },
                                    v: {
                                      type: "string",
                                    },
                                  },
                                  required: ["summary", "v"],
                                  type: "object",
                                },
                                {
                                  additionalProperties: false,
                                  properties: {
                                    rotational: {
                                      items: [
                                        {
                                          anyOf: [
                                            {
                                              additionalProperties: false,
                                              properties: {
                                                J: {
                                                  type: "number",
                                                },
                                                summary: {
                                                  type: "string",
                                                },
                                              },
                                              required: ["J", "summary"],
                                              type: "object",
                                            },
                                            {
                                              additionalProperties: false,
                                              properties: {
                                                J: {
                                                  type: "string",
                                                },
                                                summary: {
                                                  type: "string",
                                                },
                                              },
                                              required: ["J", "summary"],
                                              type: "object",
                                            },
                                          ],
                                        },
                                      ],
                                      minItems: 1,
                                      type: "array",
                                    },
                                    summary: {
                                      type: "string",
                                    },
                                    v: {
                                      type: "number",
                                    },
                                  },
                                  required: ["summary", "v"],
                                  type: "object",
                                },
                              ],
                            },
                          ],
                          minItems: 1,
                          type: "array",
                        },
                      },
                      required: ["Lambda", "S", "e", "parity", "summary"],
                      type: "object",
                    },
                  ],
                },
                type: "array",
              },
              {
                items: {
                  anyOf: [
                    {
                      additionalProperties: false,
                      properties: {
                        e: {
                          type: "string",
                        },
                        summary: {
                          type: "string",
                        },
                      },
                      required: ["e", "summary"],
                      type: "object",
                    },
                    {
                      additionalProperties: false,
                      properties: {
                        Lambda: {
                          type: "number",
                        },
                        S: {
                          type: "number",
                        },
                        e: {
                          type: "string",
                        },
                        parity: {
                          enum: ["g", "u"],
                          type: "string",
                        },
                        reflection: {
                          enum: ["+", "-"],
                          type: "string",
                        },
                        summary: {
                          type: "string",
                        },
                        vibrational: {
                          items: [
                            {
                              anyOf: [
                                {
                                  additionalProperties: false,
                                  properties: {
                                    summary: {
                                      type: "string",
                                    },
                                    v: {
                                      type: "string",
                                    },
                                  },
                                  required: ["summary", "v"],
                                  type: "object",
                                },
                                {
                                  additionalProperties: false,
                                  properties: {
                                    rotational: {
                                      items: [
                                        {
                                          anyOf: [
                                            {
                                              additionalProperties: false,
                                              properties: {
                                                J: {
                                                  type: "number",
                                                },
                                                summary: {
                                                  type: "string",
                                                },
                                              },
                                              required: ["J", "summary"],
                                              type: "object",
                                            },
                                            {
                                              additionalProperties: false,
                                              properties: {
                                                J: {
                                                  type: "string",
                                                },
                                                summary: {
                                                  type: "string",
                                                },
                                              },
                                              required: ["J", "summary"],
                                              type: "object",
                                            },
                                          ],
                                        },
                                      ],
                                      minItems: 1,
                                      type: "array",
                                    },
                                    summary: {
                                      type: "string",
                                    },
                                    v: {
                                      items: {
                                        type: "number",
                                      },
                                      type: "array",
                                    },
                                  },
                                  required: ["summary", "v"],
                                  type: "object",
                                },
                              ],
                            },
                          ],
                          minItems: 1,
                          type: "array",
                        },
                      },
                      required: ["Lambda", "S", "e", "parity", "summary"],
                      type: "object",
                    },
                  ],
                },
                type: "array",
              },
            ],
          },
          id: {
            type: "string",
          },
          particle: {
            type: "string",
          },
        },
        required: ["charge", "id", "particle"],
        type: "object",
      },
      type: "array",
    },
  },
  required: ["processes", "states"],
  type: "object",
};
